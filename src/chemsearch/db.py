import os
import logging

import pandas as pd
from rdkit.Chem import MolFromSmiles, MolFromSmarts
from rdkit.DataStructs import TanimotoSimilarity

from . import paths
from .molecule import LocalMolecule, Molecule


class MolException(Exception):
    pass


_logger = logging.getLogger(__name__)


def load_molecules(load_rdkit_mol=True):
    # drive_path = os.path.join(LOCAL_DB_PATH, 'gdrive.tsv')
    # gd = pd.read_csv(drive_path, sep='\t')
    if not os.path.exists(paths.REFERENCE_PATH):
        _logger.warning(f"Reference path not found: {paths.REFERENCE_PATH}.")
        return
    _logger.info(f"Attempting to load molecule info from {paths.REFERENCE_PATH}.")
    try:
        df = _load_reference_file_as_df()
    except pd.errors.EmptyDataError:
        return
    # categories = sorted(list(df.category.unique()))
    # category_counts = df.category.value_counts().to_dict()
    for ind, rec in df.iterrows():
        yield LocalMolecule(rec, store_mol=load_rdkit_mol)


def reload_molecules():
    global LOCAL_MOLECULES, MOLECULE_DICT
    new_molecules = list(load_molecules(load_rdkit_mol=True))
    new_id_dict = {i.inchi_key: i for i in LOCAL_MOLECULES}
    LOCAL_MOLECULES.clear()
    LOCAL_MOLECULES.extend(new_molecules)
    MOLECULE_DICT.clear()
    MOLECULE_DICT.update(new_id_dict)


def _load_reference_file_as_df():
    df = pd.read_csv(paths.REFERENCE_PATH, sep='\t', parse_dates=['mod_time'],
                     infer_datetime_format=True)
    df.sort_values(['mod_time'], inplace=True, ascending=False)
    return df


def get_substructure_matches(query_str, mols=None, is_smarts=False):
    """Get LocalMolecule matches with SMILES or SMARTS query as substructure"""
    query = MolFromSmarts(query_str) if is_smarts else MolFromSmiles(query_str)
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
