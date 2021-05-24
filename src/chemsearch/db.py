import os
import logging
from collections import Counter, defaultdict, namedtuple

import pandas as pd
from rdkit.Chem import MolFromSmiles, MolFromSmarts

from . import paths, similarity
from .molecule import LocalMolecule, Molecule, MolFileNotFoundError


class MolException(Exception):
    pass


_logger = logging.getLogger(__name__)


def iter_molecules(load_rdkit_mol=True):
    # drive_path = os.path.join(LOCAL_DB_PATH, 'gdrive.tsv')
    # gd = pd.read_csv(drive_path, sep='\t')
    if not os.path.exists(paths.REFERENCE_PATH):
        _logger.warning(f"Reference path not found: {paths.REFERENCE_PATH}.")
        return
    _logger.debug(f"Attempting to load molecule info from {paths.REFERENCE_PATH}.")
    try:
        df = _load_reference_file_as_df()
    except pd.errors.EmptyDataError:
        return
    # categories = sorted(list(df.category.unique()))
    # category_counts = df.category.value_counts().to_dict()
    for ind, rec in df.iterrows():
        yield LocalMolecule(rec, store_mol=load_rdkit_mol)


# def load_initial_db_data():
#     global LOCAL_MOLECULES, MOLECULE_DICT, DUPLICATE_TRACKER
#     try:
#         LOCAL_MOLECULES = list(iter_molecules(load_rdkit_mol=True))
#     except MolFileNotFoundError as e:
#         _logger.info(f"Aborting molecule load. {e}")
#         LOCAL_MOLECULES = []
#     MOLECULE_DICT = {i.inchi_key: i for i in LOCAL_MOLECULES}
#     DUPLICATE_TRACKER = DuplicateTracker()


def reload_molecules():
    global LOCAL_MOLECULES, MOLECULE_DICT, DUPLICATE_TRACKER
    _logger.info('Loading molecule metadata.')
    try:
        new_molecules = list(iter_molecules(load_rdkit_mol=True))
    except MolFileNotFoundError as e:
        _logger.info(f"Aborting molecule load. {e}")
        new_molecules = []
    new_id_dict = {i.inchi_key: i for i in new_molecules}
    LOCAL_MOLECULES.clear()
    LOCAL_MOLECULES.extend(new_molecules)
    MOLECULE_DICT.clear()
    MOLECULE_DICT.update(new_id_dict)
    DUPLICATE_TRACKER.update()


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
    fp_q = Molecule.get_fingerprint(query)
    results = []
    mols = LOCAL_MOLECULES if mols is None else mols
    for m in mols:
        fp = m.fingerprint_similarity_raw
        sim = similarity.calculate_similarity(fp_q, fp)
        results.append((sim, m))
    results.sort(key=lambda t: t[0], reverse=True)
    sims, molecules = zip(*results) if results else ([], [])
    return sims, molecules


def valid_mols_present():
    has_valid = False
    for mol in LOCAL_MOLECULES:
        if mol.is_valid:
            has_valid = True
            break
    return has_valid


class DuplicateTracker:

    MolId = namedtuple('MolId', ['lm_ind', 'category', 'mol_name'])

    def __init__(self, mol_list=None):
        self.has_duplicates = False
        self.duplicated_inchi_keys = set()
        self.n_duplicates = 0
        self.inchi_key_to_ids = dict()
        self.inchi_key_to_inds = dict()
        self.inchi_key_to_mols = dict()
        self.ind_to_inds = dict()
        self.ind_to_mols = dict()
        self.names_to_inds = dict()
        self.names_to_mols = dict()
        self.update(mol_list)

    def update(self, mol_list=None):
        if mol_list is None:
            mol_list = LOCAL_MOLECULES
        else:
            mol_list = list(mol_list)  # list required for indexing
        if not mol_list:
            return
        counts = Counter([i.inchi_key for i in mol_list])
        for inchi_key, count in counts.most_common():
            if count > 1:
                self.duplicated_inchi_keys.add(inchi_key)
            else:
                break
        if self.duplicated_inchi_keys:
            self.has_duplicates = True
        self.n_duplicates = len(self.duplicated_inchi_keys)
        repeats = defaultdict(set)
        for ind, i in enumerate(mol_list):
            if i.inchi_key in self.duplicated_inchi_keys:
                repeats[i.inchi_key].add(DuplicateTracker.MolId(
                    ind, i.category, i.mol_name))
        for inchi_key in self.duplicated_inchi_keys:
            id_set = repeats[inchi_key]
            inds = sorted([i.lm_ind for i in id_set])
            self.inchi_key_to_ids[inchi_key] = id_set
            self.inchi_key_to_inds[inchi_key] = inds
            self.inchi_key_to_mols[inchi_key] = [mol_list[i] for i in inds]
            for mol_id in id_set:
                matched_ids = id_set.difference((mol_id,))
                ind = mol_id.lm_ind
                matched_inds = sorted([i.lm_ind for i in matched_ids])
                matched_mols = [mol_list[i] for i in matched_inds]
                self.ind_to_inds[ind] = matched_inds
                self.ind_to_mols[ind] = matched_mols
                names = (mol_id.category, mol_id.mol_name)
                self.names_to_inds[names] = matched_inds
                self.names_to_mols[names] = matched_mols
        if self.has_duplicates:
            self.print_duplicates()

    def print_duplicates(self):
        if self.has_duplicates:
            _logger.warning("DUPLICATE INCHIKEYS FOUND!")
            for inchi_key in self.duplicated_inchi_keys:
                ids = self.inchi_key_to_ids[inchi_key]
                names = [(i.category, i.mol_name) for i in ids]
                key_str = f"{inchi_key} duplicated in: " + ', '.join([str(i) for i in names])
                _logger.info(key_str)
        else:
            _logger.info("No duplicates found.")


LOCAL_MOLECULES = []
MOLECULE_DICT = {}
DUPLICATE_TRACKER = DuplicateTracker()
