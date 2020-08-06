"""Classes for molecules (user input version and local extended info)."""
import os
import base64
from urllib.parse import urljoin

from flask import url_for
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


class Molecule:
    """Basic molecule stats object."""
    fields_stat = (
        'smiles',
        # 'smarts',
        # 'inchi',
        'inchi_key',
        # 'fingerprint_substructure',
        # 'fingerprint_similarity',
    )

    def __init__(self, mol):
        """Build stats from rdchem.Mol object."""
        self.mol = mol
        self.smiles = Chem.MolToSmiles(mol)
        # self.smarts = Chem.MolToSmarts(mol)
        # self.inchi = Chem.MolToInchi(mol)
        self.inchi_key = Chem.MolToInchiKey(mol)  # google-able. "Phenethylcyclohexane"
        # self.fingerprint_substructure = Chem.RDKFingerprint(mol).ToBase64()
        morgan_fingerprint = Molecule.get_morgan_fingerprint(mol)
        self.fingerprint_similarity_raw = morgan_fingerprint
        # morgan_base64 = base64.b64encode(morgan_fingerprint.ToBinary()).decode('utf8')
        # self.fingerprint_similarity = morgan_base64

    @classmethod
    def get_morgan_fingerprint(cls, mol):
        return Chem.AllChem.GetMorganFingerprint(mol, 2)


class LocalMolecule(Molecule):
    """Extended molecule info for local archive."""

    fields_local = (
        'mol_id',
        'mol_name',
        'mol_filename',
        'category',
        'user',
        'folder_id',
        'mod_time',
        'mol_basename',
        'dir_url',
    )

    fields_all = tuple(list(fields_local) + list(Molecule.fields_stat))

    def __init__(self, record: pd.Series, from_summary=True, store_mol=True):
        """Initialize from Google Drive molfile table or summary table record."""
        if from_summary:
            for field in self.fields_all:  # populate metadata
                self.__setattr__(field, record[field])
            self.mol_path = self._get_mol_path()
            if store_mol:
                mol = Chem.MolFromMolFile(self.mol_path)
                self.mol = mol
                self.fingerprint_similarity_raw = Molecule.get_morgan_fingerprint(mol)
            return
        # otherwise from google drive molfile table
        self.mol_path = self._get_mol_path_from_record(record)
        m = Chem.MolFromMolFile(self.mol_path)
        super().__init__(m)
        if not from_summary:  # add metadata
            self.mol_id = record['id']  # MOL file ID in
            self.mol_name = record['folder_name']
            self.mol_filename = record['name']
            self.folder_id = record['folder_id']
            self.category = record.category
            self.user = record.lastModifyingUser
            self.mod_time = record.modifiedTime
            self.mol_basename = os.path.splitext(self.mol_filename)[0]
            self.dir_url = f"https://drive.google.com/drive/u/0/folders/{self.folder_id}"

    @property
    def url_svg(self):
        data_dir = url_for('static', filename='data')
        local_dir = '/'.join([data_dir, self.category, self.mol_name])
        svg_path = '/'.join([local_dir, 'ref_image.svg'])
        svg_url = urljoin(data_dir, svg_path)
        return svg_url

    def _get_mol_path(self):
        archive_dir = os.environ.get('LOCAL_DB_PATH', 'local_db')
        mol_dir = os.path.join(archive_dir, self.category, self.mol_name)
        mol_path = os.path.join(mol_dir, self.mol_filename)
        return mol_path

    @staticmethod
    def _get_mol_path_from_record(record=None):
        """Get mol path from Drive table record."""
        archive_dir = os.environ.get('LOCAL_DB_PATH', 'local_db')
        mol_dir = os.path.join(archive_dir, record.category, record.folder_name)
        mol_path = os.path.join(mol_dir, record['name'])
        return mol_path

    def __repr__(self):
        return f"<Molecule {self.mol_id}>"
