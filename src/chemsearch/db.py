import os
import logging

import pandas as pd
from rdkit.Chem import MolFromSmiles
from rdkit.DataStructs import TanimotoSimilarity

from .molecule import LocalMolecule, Molecule
from .admin import REFERENCE_PATH


class MolException(Exception):
    pass


_logger = logging.getLogger(__name__)


def load_molecules(load_rdkit_mol=True):
    # drive_path = os.path.join(LOCAL_DB_PATH, 'gdrive.tsv')
    # gd = pd.read_csv(drive_path, sep='\t')
    if not os.path.exists(REFERENCE_PATH):
        _logger.warning(f"Reference path not found: {REFERENCE_PATH}.")
        return
    _logger.info(f"Attempting to load molecule info from {REFERENCE_PATH}.")
    try:
        df = pd.read_csv(REFERENCE_PATH, sep='\t', parse_dates=['mod_time'],
                         infer_datetime_format=True)
        df.sort_values(['mod_time'], inplace=True, ascending=False)
    except pd.errors.EmptyDataError:
        return
    # categories = sorted(list(df.category.unique()))
    # category_counts = df.category.value_counts().to_dict()
    for ind, rec in df.iterrows():
        yield LocalMolecule(rec, store_mol=load_rdkit_mol)


def reload_molecules():
    global LOCAL_MOLECULES
    LOCAL_MOLECULES.clear()
    LOCAL_MOLECULES.extend(load_molecules(load_rdkit_mol=True))


def get_single_molecule():
    gen = load_molecules()
    mol = next(gen)
    return mol


def get_substructure_matches(smiles, mols=None):
    """Get LocalMolecule matches with smiles as substructure"""
    query = MolFromSmiles(smiles)
    if query is None:
        raise MolException("Error building molecule from input.")
    mols = LOCAL_MOLECULES if mols is None else mols
    matches = [i for i in mols if i.mol.HasSubstructMatch(query)]
    return matches


def get_sim_matches(smiles, mols=None):
    query = MolFromSmiles(smiles)
    if query is None:
        raise MolException("Error building molecule from input.")
    fp_q = Molecule.get_morgan_fingerprint(query)
    results = []
    mols = LOCAL_MOLECULES if mols is None else mols
    for m in mols:
        fp = m.fingerprint_similarity_raw
        sim = TanimotoSimilarity(fp_q, fp)
        results.append((sim, m))
    results.sort(key=lambda t: t[0], reverse=True)
    sims, molecules = zip(*results) if results else ([], [])
    return sims, molecules


LOCAL_MOLECULES = list(load_molecules(load_rdkit_mol=True))
MOLECULE_DICT = {i.inchi_key: i for i in LOCAL_MOLECULES}
