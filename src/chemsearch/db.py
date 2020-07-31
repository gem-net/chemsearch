import os

import pandas as pd
from rdkit.Chem import MolFromSmiles
from rdkit.DataStructs import TanimotoSimilarity

from .models import LocalMolecule, Molecule


LOCAL_DB_PATH = os.environ.get('LOCAL_DB_PATH')


class MolException(Exception):
    pass


def get_molecules(load_rdkit_mol=True):
    # drive_path = os.path.join(LOCAL_DB_PATH, 'gdrive.tsv')
    # gd = pd.read_csv(drive_path, sep='\t')
    summary_path = os.path.join(LOCAL_DB_PATH, 'summary.tsv')
    df = pd.read_csv(summary_path, sep='\t')
    # categories = sorted(list(df.category.unique()))
    # category_counts = df.category.value_counts().to_dict()
    for ind, rec in df.iterrows():
        yield LocalMolecule(rec, store_mol=load_rdkit_mol)


def get_single_molecule():
    gen = get_molecules()
    mol = next(gen)
    return mol


def get_substructure_matches(smiles):
    """Get LocalMolecule matches with smiles as substructure"""
    query = MolFromSmiles(smiles)
    if query is None:
        raise MolException("Error building molecule from input.")
    matches = [i for i in get_molecules() if i.mol.HasSubstructMatch(query)]
    return matches


def get_sim_matches(smiles):
    query = MolFromSmiles(smiles)
    if query is None:
        raise MolException("Error building molecule from input.")
    fp_q = Molecule.get_morgan_fingerprint(query)
    results = []
    for m in get_molecules(load_rdkit_mol=True):
        fp = m.fingerprint_similarity_raw
        sim = TanimotoSimilarity(fp_q, fp)
        results.append((sim, m))

    results.sort(key=lambda t: t[0], reverse=True)
    sims, molecules = zip(*results)
    return sims, molecules
