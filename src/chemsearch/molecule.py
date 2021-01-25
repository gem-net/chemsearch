"""Classes for molecules (user input version and local extended info)."""
import os
# import base64
import logging
from urllib.parse import urljoin

from flask import url_for
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from . import paths


_logger = logging.getLogger(__name__)


class Molecule:
    """Basic molecule stats object."""
    fields_stat = (
        'is_valid',
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
        if mol is not None:
            self.is_valid = True
            self.smiles = Chem.MolToSmiles(mol)
            # self.smarts = Chem.MolToSmarts(mol)
            # self.inchi = Chem.MolToInchi(mol)
            self.inchi_key = Chem.MolToInchiKey(mol)  # google-able. "Phenethylcyclohexane"
            # self.fingerprint_substructure = Chem.RDKFingerprint(mol).ToBase64()
            self.fingerprint_similarity_raw = Molecule.get_morgan_fingerprint(mol)
            # morgan_base64 = base64.b64encode(morgan_fingerprint.ToBinary()).decode('utf8')
            # self.fingerprint_similarity = morgan_base64
        else:
            self.is_valid = False
            self.smiles = None
            self.inchi_key = None
            self.fingerprint_similarity_raw = None

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
        """Initialize from scanned molfile table or summary table record."""
        if from_summary:  # all metadata fields in record, via REFERENCE_PATH
            for field in self.fields_all:  # populate metadata
                self.__setattr__(field, record[field])
            self.mol_path = self._get_mol_path()
            if store_mol:
                self.mol = LocalMolecule._load_structure_from_mol(self.mol_path)
                if self.mol is not None:
                    self.fingerprint_similarity_raw = Molecule.get_morgan_fingerprint(self.mol)
                    self.mol.SetProp("_Name", self.mol_name)
                else:
                    self.fingerprint_similarity_raw = None
            return
        # otherwise record is from SCAN_RESULTS
        self.mol_path = self._get_mol_path_from_scan_record(record)
        m = LocalMolecule._load_structure_from_mol(self.mol_path)
        super().__init__(m)
        self.mol_id = record['id']
        self.mol_name = record['folder_name']
        self.mol_filename = record['name']
        self.folder_id = record['folder_id']
        self.category = record['category']
        self.user = record['lastModifyingUser']
        self.mod_time = record['modifiedTime']
        self.mol_basename = os.path.splitext(self.mol_filename)[0]
        self.dir_url = f"https://drive.google.com/drive/u/0/folders/{self.folder_id}"
        if self.is_valid:
            self.mol.SetProp("_Name", self.mol_name)

    @property
    def url_svg(self):
        data_dir = url_for('static', filename='data')
        local_dir = '/'.join([data_dir, self.category, self.mol_name])
        svg_path = '/'.join([local_dir, self._get_svg_basename()])
        svg_url = urljoin(data_dir, svg_path)
        return svg_url

    def _get_mol_path(self):
        mol_dir = self.local_mol_dir
        mol_path = os.path.join(mol_dir, self.mol_filename)
        return mol_path

    @property
    def local_mol_dir(self):
        return os.path.join(paths.ARCHIVE_DIR, self.category, self.mol_name)

    def _get_svg_basename(self):
        return f'ref_{self.inchi_key}.svg'

    @property
    def svg_path(self):
        svg_basename = self._get_svg_basename()
        return os.path.join(self.local_mol_dir, svg_basename)

    @staticmethod
    def _get_mol_path_from_scan_record(record=None):
        """Get mol path from SCAN_RESULTS record (Drive or local)."""
        mol_dir = os.path.join(paths.ARCHIVE_DIR, record.category,
                               record.folder_name)
        mol_path = os.path.join(mol_dir, record['name'])
        return mol_path

    @staticmethod
    def _load_structure_from_mol(mol_path):
        if not os.path.exists(mol_path):
            raise MolFileNotFoundError(f"{mol_path} not found.")
        mol = Chem.MolFromMolFile(mol_path)
        return mol

    def __repr__(self):
        return f"<Molecule {self.mol_id}>"


class MolFileNotFoundError(Exception):
    pass
