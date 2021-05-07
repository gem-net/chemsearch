"""Assembles metadata (SMILES, inchi, fingerprints) and reference images for
local archive."""

import os
import glob
import logging
import datetime
import pathlib
import hashlib
import shutil
import tarfile
from importlib import resources
from collections import OrderedDict

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole as ipc
from IPython.display import display
from dotenv import load_dotenv, find_dotenv

from .molecule import LocalMolecule
from .plot import save_svg_if_not_present
from . import drive, paths

_logger = logging.getLogger(__name__)


def create_support_dirs_extract_resources():
    os.makedirs(paths.DATA_ROOT, exist_ok=True)
    os.makedirs(paths.DEMO_DIR, exist_ok=True)
    os.makedirs(paths.CONFIG_DIR, exist_ok=True)
    # Extract demo_db into user_data_dir
    with resources.path('chemsearch', 'demo_db.tar.gz') as gz_path:
        with tarfile.open(gz_path, "r:gz") as tar:
            _logger.info(f"Extracting demo data to {paths.DATA_ROOT}.")
            tar.extractall(path=paths.DATA_ROOT)
    # Extract config files
    added_config = []
    for config_basename in ['demo.env',
                            'demo.custom_queries.yaml',
                            'external_dbs.yaml']:
        new_path = paths.CONFIG_DIR.joinpath(config_basename)
        if not new_path.exists():
            added_config.append(config_basename)
            with resources.path('chemsearch.config', config_basename) as path:
                shutil.copy(path, new_path)
    if added_config:
        _logger.info(f"Copied config files to {paths.CONFIG_DIR}: {added_config}.")
    else:
        _logger.info(f"Config files found in {paths.CONFIG_DIR}.")


def run_full_scan_and_rebuild_from_drive_no_app():
    meta = drive.Meta().build()
    reload_env()
    drive.create_local_archive(meta.molfiles, local_root=paths.ARCHIVE_DIR,
                               files_resource=meta.files_resource,
                               scan_path=paths.SCAN_RESULTS_PATH)
    df = assemble_archive_metadata(paths.ARCHIVE_DIR)
    return df


def scan_local_archive():
    """Builds local.tsv for local MOL files, equiv of gdrive.tsv Drive info."""
    root = pathlib.Path(paths.ARCHIVE_DIR)
    records = []
    root_dir, test_categs, _ = next(os.walk(root, topdown=True))
    for test_categ in test_categs:
        categ_path = root.joinpath(test_categ)
        _, test_mol_names, _ = next(os.walk(categ_path))
        for test_mol_name in test_mol_names:
            moldir_path = categ_path.joinpath(test_mol_name)
            # Valid if moldir has MOL file
            _, _, files = next(os.walk(moldir_path))
            mol_files = [i for i in files if i.endswith('.mol')]
            if len(mol_files) != 1:
                _logger.warning(f"Skipping {test_categ} {test_mol_name}: {len(mol_files)} mol files found.")
                continue
            mol_file = mol_files[0]
            mol_path = moldir_path.joinpath(mol_file)
            md5 = _get_md5(mol_path)
            mtime = _get_mtime(mol_path)
            mtime_dir = _get_mtime(moldir_path)
            record = OrderedDict({
                'id': md5,
                'name': mol_file,
                'modifiedTime': mtime,
                'lastModifyingUser': None,
                'category': test_categ,
                'folder_id': None,
                'folder_name': test_mol_name,
                'folder_modified': mtime_dir,
                'folder_user': None,
            })
            records.append(record)
    if records:
        df = pd.DataFrame(data=records)
    else:
        df = pd.DataFrame(
            columns=['id', 'name', 'modifiedTime', 'lastModifyingUser', 'category',
                     'folder_id', 'folder_name', 'folder_modified', 'folder_user'])
    df.to_csv(paths.SCAN_RESULTS_PATH, sep='\t', index=False)
    return df


def assemble_archive_metadata(archive_dir=None, use_drive=False):
    """Assemble molecule associated data, inc images, fingerprints.

    Export summary table (summary.tsv) to archive directory root.
    """
    paths.update_paths(use_drive=use_drive)
    mols = pd.read_csv(paths.SCAN_RESULTS_PATH, sep='\t',
                       parse_dates=['modifiedTime', 'folder_modified'],
                       infer_datetime_format=True)
    mol_info = []
    for ind, mol in mols.iterrows():
        info = OrderedDict()
        mol_dir = os.path.join(archive_dir, mol.category, mol.folder_name)
        _logger.info(f"Processing directory: {mol_dir}")
        m = LocalMolecule(mol, from_summary=False)
        if m.is_valid:
            save_svg_if_not_present(m.mol, m.svg_path)
        else:
            _logger.warning(f"Invalid MOL for {m.mol_name}.")
        for field in m.fields_all:
            info[field] = getattr(m, field)
        mol_info.append(info)
    summary = pd.DataFrame.from_records(mol_info)
    summary.to_csv(paths.REFERENCE_PATH, sep='\t', index=False)
    _logger.info(f"Reference images written to molecule directories.")
    _logger.info(f"Molecule identifiers and fingerprints written to {paths.REFERENCE_PATH}.")
    return summary


def _get_mtime(path):
    stats = os.stat(path)
    return datetime.datetime.fromtimestamp(stats.st_mtime)


def _get_md5(file_path):
    h = hashlib.md5()
    with open(file_path, 'rb') as infile:
        for line in infile:
            h.update(line)
    return h.hexdigest()


def save_mol_image_ipc(mol, im_path):
    ipc.ipython_useSVG = True  # True
    im_format = os.path.splitext(im_path)[1][1:]
    fn_mode = {'png': (ipc._toPNG, 'wb'),
               'svg': (ipc._toSVG, 'w')}
    fn, mode = fn_mode[im_format]
    imdata = fn(mol)
    with open(im_path, mode) as out:
        out.write(imdata)


def demo_gather_metadata_stage_dir(archive_dir=None):
    """Gather metadata for molecules in subdirectories of staging area."""

    mol_dirs = []
    categ_dirs = [i.path for i in os.scandir(archive_dir) if i.is_dir()]
    for categ_dir in categ_dirs:
        mol_dirs.extend([i.path for i in os.scandir(categ_dir) if i.is_dir()])

    for mol_dir in mol_dirs:
        mol_paths = glob.glob(os.path.join(mol_dir, '*.mol'))

        if len(mol_paths) != 1:
            _logger.warning(f"Skipping {mol_dir}: found {len(mol_paths)} MOL files. Need 1.")
            continue

        mol_path = mol_paths[0]
        mol_dir, mol_basename = os.path.split(mol_path)

        _logger.info(f"Organizing molecule directory: {mol_dir}")

        m = Chem.MolFromMolFile(mol_path)
        display(m)

        smiles = Chem.MolToSmiles(m)
        _logger.info(f"  smiles: {smiles}")

        smarts = Chem.MolToSmarts(m)
        _logger.info(f"  smarts: {smarts}")

        inchi = Chem.MolToInchi(m)
        _logger.info(f"  inchi: {inchi}")

        inchi_key = Chem.MolToInchiKey(m)  # google-able. "Phenethylcyclohexane"
        _logger.info(f"  inchi_key: {inchi_key}")

        # m2 = Chem.MolFromSmiles(smarts)

        Draw.MolToFile(m, os.path.join(mol_dir, 'ref_image.png'))
        Draw.MolToFile(m, os.path.join(mol_dir, 'ref_image.svg'))
        _logger.info("  Reference images written to png and svg.\n")


def reload_env():
    _logger.info("Loading .env from %s", paths.ENV_PATH)
    load_dotenv(paths.ENV_PATH, override=True, verbose=True)
