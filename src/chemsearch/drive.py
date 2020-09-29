"""Queries Drive API to identify categories and MOL files. Builds local archive."""
import io
import os
import time
import logging

import pandas as pd
from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload, HttpError
from bs4 import BeautifulSoup

from .helpers import parse_timestamp_str


_logger = logging.getLogger(__name__)
GDOC_MIMETYPE = 'application/vnd.google-apps.document'
TXT_MIMETYPE = 'text/plain'

MIME_MAP = {
    GDOC_MIMETYPE: ('text/html', '.html'),  # ('application/pdf', '.pdf'),  # Google Docs
    'application/vnd.google-apps.spreadsheet':
        ('application/vnd.openxmlformats-officedocument.spreadsheetml.sheet', '.xlsx'),  # Google Sheets
    'application/vnd.google-apps.presentation': ('application/pdf', '.pdf'),  # Google Slides
}


class Meta:
    """Hold file metadata for molecules in Shared Drive."""
    def __init__(self):
        """Create empty metadata object."""
        self.category_dict = dict()
        self.folders = pd.DataFrame()
        self.molfiles = pd.DataFrame()
        self._files_resource = None

    def build(self):
        """Gather all metadata using Drive API."""
        self._load_files_resource()
        categories = get_category_ids(files_resource=self.files_resource)
        _logger.info(f"Identified {len(categories)} categories.")
        folders = get_mol_folders(categories, files_resource=self.files_resource)
        _logger.info(f"Identified {len(folders)} molecule directories.")
        mols = get_mol_files(folders, files_resource=self.files_resource)
        _logger.info(f"Identified {len(mols)} MOL files.")
        self.category_dict = categories
        self.folders = folders
        self.molfiles = mols
        return self

    @property
    def latest_mol_time(self):
        if not len(self.molfiles):
            return None
        last_mod = self.molfiles.modifiedTime.sort_values(ascending=False).iloc[0]
        return int(last_mod.timestamp())

    @property
    def files_resource(self):
        if self._files_resource is None:
            self._load_files_resource()
        return self._files_resource

    def _load_files_resource(self):
        resource = get_files_service().files()
        self._files_resource = resource


def create_local_archive(mols, local_root=None, files_resource=None,
                         scan_path=None):
    """Download all MOL files to local directory (unless they already exist).

    Local archive has category folders, with nested molecule directories.

    Args:
        mols (pd.DataFrame): MOL table, output from get_mol_files.
        local_root (str): path to directory in which to store files.
        files_resource: files API resource
        scan_path (str): TSV output path for Drive MOL file summary.
    """
    if not local_root:
        local_root = os.path.join(os.getcwd(), 'local_db')
        os.makedirs(local_root, exist_ok=True)
        _logger.info(f"Using {local_root} for local archive.")
    else:
        local_root = os.path.abspath(local_root)
    to_download = []
    to_skip = []
    for ind, mol in mols.iterrows():
        basename = mol['name']
        mol_id = mol['id']
        # url = mol.webContentLink
        category = mol.category
        folder_name = mol.folder_name
        folder_path = os.path.join(local_root, category, folder_name)
        mol_path = os.path.join(folder_path, basename)
        if os.path.exists(mol_path):
            to_skip.append(mol_id)
        else:
            os.makedirs(folder_path, exist_ok=True)
            to_download.append((mol_id, mol_path))
    for mol_id, save_path in to_download:
        download_file(mol_id, save_path, files_resource=files_resource)
    if len(to_skip):
        _logger.info(f"Skipped {len(to_skip)} MOL files already in archive.")
    if not scan_path:
        scan_path = os.path.join(local_root, 'scan_gdrive.tsv')  # output table
    mols.to_csv(scan_path, sep='\t', index=False)


def get_mol_files(folders, files_resource=None):
    """Find MOL files and join to folder info.

    Args:
        folders (pd.DataFrame): output from get_mol_folders.
        files_resource: files API resource.
    """
    mols_full = run_query("mimeType = 'chemical/x-mdl-molfile'  and trashed = false",
                          files_resource=files_resource)  # type: pd.DataFrame
    if mols_full is None:
        return pd.DataFrame(
            columns=['id', 'name', 'modifiedTime', 'lastModifyingUser', 'category',
                     'folder_id', 'folder_name', 'folder_modified', 'folder_user'])
    mols = mols_full[['id', 'name', 'modifiedTime', 'lastModifyingUser',
                      'parents', 'webContentLink']].copy()
    mols['folder_id'] = mols['parents'].transform(lambda v: v[0])
    mols.drop('parents', axis=1, inplace=True)
    mol_info = mols.merge(folders, how='left', on='folder_id')
    # ignore mol files that aren't in top molecule folders
    n_ignore = mol_info.category.isnull().sum()
    if n_ignore:
        _logger.info(f"Skipping {n_ignore} molecules in subdirectories.")
        mol_info.dropna(axis=0, subset=['category'], inplace=True)
    return mol_info


def get_category_ids(files_resource=None):
    """Get all molecule categories and their associated folder IDs.

    Returns:
        dictionary of category name -> folder ID.
    """
    compounds_root = os.environ.get('SHARED_DRIVE_ID')
    query = f"{compounds_root!r} in parents and trashed = false and mimeType = 'application/vnd.google-apps.folder'"
    files = run_query(query, as_df=False, files_resource=files_resource)
    categories = {file['name']: file['id'] for file in files}
    return categories


def get_mol_folders(category_dict, files_resource=None):
    """Get table of folder metadata for all molecule folders in all categories."""
    t1 = time.perf_counter()
    df_list = []
    for category in category_dict:
        temp = scan_folder(category_dict[category], files_resource=files_resource)
        if temp is not None:
            temp.insert(0, 'category', category)
            df_list.append(temp)
    t2 = time.perf_counter()
    _logger.info(f"Folder lookup complete in {t2 - t1:.1f} seconds.")
    if not df_list:
        return pd.DataFrame()
    folders = pd.concat(df_list, axis=0, ignore_index=True)
    folder_rename = {
        'id': 'folder_id',
        'name': 'folder_name',
        'modifiedTime': 'folder_modified',
        'lastModifyingUser': 'folder_user',
    }
    folders = folders.rename(columns=folder_rename)[
        ['category'] + [i for i in folder_rename.values()]]
    return folders


def scan_folder(folder_id, file_fields=None, files_resource=None):
    """Get table of folder contents, excluding trashed items."""
    query = "{!r} in parents and trashed = false".format(folder_id)
    df = run_query(query, file_fields=file_fields,
                   files_resource=files_resource)
    return df


def run_query(query, page_size=1000, as_df=True, file_fields=None,
              files_resource=None):
    """Search Google Drive using query."""
    # folder_id = os.environ.get('COMPOUNDS_DIR_ID')
    team_drive_id = os.environ.get('SHARED_DRIVE_ID')
    if not files_resource:
        files_resource = get_files_service().files()
    if not file_fields:
        file_cols = ['name', 'id', 'kind', 'parents', 'mimeType', 'createdTime',
                     'modifiedTime', 'trashed', 'explicitlyTrashed',
                     'lastModifyingUser/displayName', 'webContentLink',
                     'iconLink', 'webViewLink']
        file_fields = f"nextPageToken, files({', '.join(file_cols)})"

    page_token = None
    file_list = []
    while True:
        res = files_resource.list(q=query,
                                  corpora='drive', spaces='drive',
                                  includeItemsFromAllDrives=True,
                                  orderBy='modifiedTime desc',
                                  pageSize=page_size,
                                  supportsAllDrives=True,
                                  driveId=team_drive_id,
                                  fields=file_fields,
                                  pageToken=page_token).execute()
        file_list.extend(res.get('files', []))
        page_token = res.get('nextPageToken', None)
        if page_token is None:
            break
    if not as_df:
        return file_list
    files = pd.DataFrame.from_records(file_list)  # type: pd.DataFrame
    if len(files):
        files['lastModifyingUser'] = files['lastModifyingUser'].transform(lambda v: v['displayName'])
        files['createdTime'] = files['createdTime'].apply(parse_timestamp_str)
        files['modifiedTime'] = files['modifiedTime'].apply(parse_timestamp_str)
        files['is_folder'] = files.mimeType.apply(lambda v: v.endswith('folder'))
        return files


def get_files_service():
    """Get dictionary of {service_name: service_handle}."""
    service_account_file = os.getenv('SERVICE_ACCOUNT_FILE')
    scopes = ['https://www.googleapis.com/auth/drive.readonly',
              'https://www.googleapis.com/auth/drive.metadata.readonly']

    credentials = service_account.Credentials.from_service_account_file(
        service_account_file, scopes=scopes)
    delegated_credentials = credentials.with_subject(
        os.environ.get('CREDENTIALS_AS_USER'))

    files_service = build('drive', 'v3', credentials=delegated_credentials,
                          cache_discovery=False)
    return files_service


def download_file(file_id, save_path, files_resource=None):
    """Get binary data for file."""
    if not files_resource:
        files_resource = get_files_service().files()
    request = files_resource.get_media(fileId=file_id)
    # fh = io.BytesIO()
    _logger.info(f"Saving MOL file to {save_path}...")
    with open(save_path, 'wb') as fh:
        downloader = MediaIoBaseDownload(fh, request)
        done = False
        while done is False:
            status, done = downloader.next_chunk()


def get_file_listing_and_custom_info(mol, files_resource=None):
    """Get molecule top-level folder listing and custom file info.

    Returns:
        df (pd.DataFrame, None): metadata for files in top-level molecule folder.
        rec (pd.Series, None): df row corresponding to custom file.
        content (html, None): content of custom file, in html.
    """
    df = scan_folder(mol.folder_id, files_resource=files_resource)
    df.sort_values('modifiedTime', ascending=False)
    if type(df) is not pd.DataFrame:
        _logger.info("Folder lookup failed for %s.", mol)
        return None, None, None
    rec, content = _get_custom_record_and_content_from_folder_df(df)
    if rec is None:
        _logger.info("No custom file found.")
        return df, None, None
    return df, rec, content


def _get_custom_record_and_content_from_folder_df(df):
    content, content_type = None, None
    rec = _get_custom_file_record_if_exists(df)  # pd.Series or None
    if type(rec) is not pd.Series:
        return None, None
    file_id, mimetype = rec['id'], rec['mimeType']
    try:
        if mimetype == GDOC_MIMETYPE:
            content = _get_doc_html_body_string(file_id)
        elif mimetype == TXT_MIMETYPE:
            content = _get_txt_file_as_string(file_id)
            content = f"<pre>{content}</pre>"
    except HttpError:
        _logger.info(f"Failed file lookup for record %s", rec.to_dict())
        return None, None
    return rec, content


def _get_custom_file_record_if_exists(df):
    is_gdoc = df['mimeType'] == GDOC_MIMETYPE
    name_is_custom = df['name'].apply(lambda v: os.path.splitext(v)[0].lower() == 'custom')

    gdoc_match = is_gdoc & name_is_custom
    if gdoc_match.sum() == 1:
        rec = df[gdoc_match].iloc[0]
        return rec
    is_txt = df['mimeType'] == TXT_MIMETYPE
    txt_match = is_txt & name_is_custom
    if txt_match.sum() == 1:
        rec = df[txt_match].iloc[0]
        return rec


def _get_file_bytes(file_id, mime_orig, title=None, files_resource=None):
    """Download Drive file, converting format if necessary.

    Args:
        file_id (str): Drive file id.
        title (str): Drive file title.
        mime_orig (str): Mime type of Drive file.
    Returns:
        fh (io.BytesIO): file stream
        filename: file and extension of output file
        mime_out: output mime type
    """
    if not files_resource:
        files_resource = get_files_service().files()
    if title is None:
        title = file_id
    if mime_orig in MIME_MAP:
        mime_out, extension = MIME_MAP[mime_orig]
        request = files_resource.export_media(fileId=file_id, mimeType=mime_out)
        filename = ''.join([title, extension])
    else:  # direct download
        request = files_resource.get_media(fileId=file_id)
        filename = title
        mime_out = mime_orig

    fh = io.BytesIO()
    downloader = MediaIoBaseDownload(fh, request)
    done = False
    while done is False:
        status, done = downloader.next_chunk()
        # print("Download %d%%." % int(status.progress() * 100))
    _logger.debug('Downloaded {}'.format(title))
    return fh, filename, mime_out


def _get_doc_html_body_string(file_id):
    fh, _, _ = _get_file_bytes(file_id, mime_orig=GDOC_MIMETYPE)
    content = fh.getvalue().decode()
    body = BeautifulSoup(content, features="lxml").body
    content = '\n'.join([str(i) for i in body.contents])
    return content


def _get_txt_file_as_string(file_id):
    fh, _, _ = _get_file_bytes(file_id, mime_orig=TXT_MIMETYPE)
    content = fh.getvalue().decode().strip()
    return content
