import io
import os
import pathlib
import logging
import datetime
from collections import OrderedDict
from itertools import product

import markdown
import pandas as pd

from ..helpers import clean_html

_logger = logging.getLogger(__name__)


def get_file_listing_and_custom_info(mol):
    """Get molecule top-level folder listing and custom file info.

    Returns:
        df (pd.DataFrame, None): metadata for files in top-level molecule folder.
        rec (pd.Series, None): df row corresponding to custom file.
        content (html, None): content of custom file, in html.
    """
    dir_path = pathlib.Path(mol.local_mol_dir)
    df = _get_listing_for_local_dir(dir_path)
    if type(df) is not pd.DataFrame:
        _logger.info("Folder lookup failed for %s.", mol)
        return None, None, None
    rec = _get_custom_file_record_if_exists(df)  # pd.Series or None
    if type(rec) is not pd.Series:
        return df, None, None
    custom_path = dir_path.joinpath(rec['name'])
    is_md = os.path.splitext(custom_path)[1].lower() == '.md'
    if is_md:
        content = _get_clean_html_from_md_file(custom_path)
    else:
        with open(custom_path, 'r') as infile:
            content = infile.read()
        content = f"<pre>{content}</pre>"
    return df, rec, content


def _get_clean_html_from_md_file(md_path):
    md_bytes = io.BytesIO()
    markdown.markdownFromFile(input=str(md_path), output_format='html',
                              output=md_bytes, tab_length=4,
                              extensions=['sane_lists', 'tables']
                              )
    html = md_bytes.getvalue().decode('utf8')
    html_clean = clean_html(html)
    return html_clean


def _get_mtime(path):
    stats = os.stat(path)
    return datetime.datetime.fromtimestamp(stats.st_mtime)


def _get_listing_for_local_dir(dir_path):
    dir_path = pathlib.Path(dir_path)
    _, subdirs, filenames = next(os.walk(dir_path))
    skip_files = frozenset(['.DS_Store', ])
    filenames = [i for i in filenames if i not in skip_files]
    types_names = list(product(['folder'], subdirs)) + list(product(['file'], filenames))
    records = []
    for kind, name in types_names:
        path = dir_path.joinpath(name)
        mtime = _get_mtime(path)
        record = OrderedDict({
            'name': name,
            'kind': kind,
            'modifiedTime': mtime,
        })
        records.append(record)
    df = pd.DataFrame(data=records) if records else None
    df.sort_values('modifiedTime', ascending=False, inplace=True)
    return df


def _get_custom_file_record_if_exists(df):
    """Priority: custom.md,  *.md x1, custom.* x1, custom.txt, .txt x1, None."""
    is_file = df['kind'] == 'file'
    name_start = df['name'].apply(lambda v: os.path.splitext(v)[0].lower())
    ext = df['name'].apply(lambda v: os.path.splitext(v)[1].lower())
    name_is_custom = name_start == 'custom'
    is_md = ext == '.md'
    is_txt = ext == '.txt'
    tests = [
        name_is_custom & is_md & is_file,
        is_md & is_file,
        name_is_custom & is_file,
        name_is_custom & is_txt & is_file,
        is_txt & is_file,
    ]
    for matches in tests:
        n_matches = matches.sum()
        if n_matches == 1:
            return df[matches].iloc[0]
