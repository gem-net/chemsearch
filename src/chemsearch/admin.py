"""Assembles metadata (SMILES, inchi, fingerprints) and reference images for
local archive."""

import os
import glob
import logging
from collections import OrderedDict

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole as ipc
from IPython.display import display
from dotenv import load_dotenv, find_dotenv

from .models import LocalMolecule
from .plot import plot_mol, save_images

_logger = logging.getLogger(__name__)


def assemble_archive_metadata(archive_dir=None):
    """Assemble molecule associated data, inc images, fingerprints.

    Export summary table (summary.tsv) to archive directory root.
    """
    if archive_dir is None:
        os.environ.get('LOCAL_DB_PATH', 'local_db')
    drive_data_path = os.path.join(archive_dir, 'gdrive.tsv')
    mols = pd.read_csv(drive_data_path, sep='\t')
    mol_info = []
    for ind, mol in mols.iterrows():
        info = OrderedDict()
        mol_dir = os.path.join(archive_dir, mol.category, mol.folder_name)
        _logger.info(f"Processing directory: {mol_dir}")
        m = LocalMolecule(mol, from_summary=False)
        save_images(m.mol, mol_dir)
        for field in m.fields_all:
            info[field] = getattr(m, field)
        mol_info.append(info)
    summary = pd.DataFrame.from_records(mol_info)
    out_path = os.path.join(archive_dir, 'summary.tsv')
    summary.to_csv(out_path, sep='\t', index=False)
    _logger.info(f"Reference images written to molecule directories.")
    _logger.info(f"Molecule identifiers and fingerprints written to {out_path}.")
    return summary


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
    env_path = os.getenv('ENV_NAME', find_dotenv(usecwd=True))
    _logger.info("Loading .env from %s", env_path)
    load_dotenv(env_path, override=True, verbose=True)
