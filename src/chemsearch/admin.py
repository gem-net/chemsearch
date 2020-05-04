"""Assembles metadata (SMILES, inchi, fingerprints) and reference images for
local archive."""

import os
import glob
import base64
import logging
from collections import OrderedDict

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from dotenv import load_dotenv, find_dotenv
from IPython.display import display

TEST_DIR = '/Users/sgg/Desktop/test_molecules'

_logger = logging.getLogger(__name__)


class MolException(Exception):
    pass


def assemble_archive_metadata(archive_dir=None):
    """Assemble molecule associated data, inc images, fingerprints.

    Export summary table (summary.tsv) to archive directory root.
    """

    drive_data_path = os.path.join(archive_dir, 'gdrive.tsv')
    mols = pd.read_csv(drive_data_path, sep='\t')
    mol_info = []
    for ind, mol in mols.iterrows():
        info = OrderedDict()
        mol_dir = os.path.join(archive_dir, mol.category, mol.folder_name)
        mol_path = os.path.join(mol_dir, mol['name'])
        _logger.info(f"Processing directory: {mol_dir}")
        m = Chem.MolFromMolFile(mol_path)
        info['molecule_id'] = mol['id']  # MOL file ID in
        info['category'] = mol.category
        info['user'] = mol.lastModifyingUser
        info['mod_time'] = mol.modifiedTime
        info['smiles'] = Chem.MolToSmiles(m)
        info['smarts'] = Chem.MolToSmarts(m)
        info['inchi'] = Chem.MolToInchi(m)
        info['inchi_key'] = Chem.MolToInchiKey(m)  # google-able. "Phenethylcyclohexane"
        # SAVE REFERENCE IMAGES
        png_path = os.path.join(mol_dir, 'ref_image.png')
        svg_path = os.path.join(mol_dir, 'ref_image.svg')
        if not os.path.exists(png_path) and not os.path.exists(svg_path):
            Chem.Draw.MolToFile(m, png_path)
            Chem.Draw.MolToFile(m, svg_path)
        info['fingerprint_substructure'] = Chem.RDKFingerprint(m).ToBase64()
        morgan_fingerprint = Chem.AllChem.GetMorganFingerprint(m, 2)
        morgan_base64 = base64.b64encode(morgan_fingerprint.ToBinary()).decode('utf8')
        info['fingerprint_similarity'] = morgan_base64
        mol_info.append(info)
    summary = pd.DataFrame.from_records(mol_info)
    out_path = os.path.join(archive_dir, 'summary.tsv')
    summary.to_csv(out_path, sep='\t', index=False)
    _logger.info(f"Reference images written to molecule directories.")
    _logger.info(f"Molecule identifiers and fingerprints written to {out_path}.")
    return summary


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
