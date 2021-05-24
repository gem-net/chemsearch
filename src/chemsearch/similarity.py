from rdkit import DataStructs
from rdkit.Chem.MolDb import FingerprintUtils as fputils


def get_gobbi_fp(mol):
    from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
    factory = Gobbi_Pharm2D.factory
    fp = Generate.Gen2DFingerprint(mol, factory)
    return fp


def get_maccs_keys(mol):
    from rdkit.Chem import MACCSkeys
    fp = MACCSkeys.GenMACCSKeys(mol)
    return fp


fp_fn_dict = {
    'Morgan': fputils.BuildMorganFP,
    'RDK': fputils.BuildRDKitFP,
    'AtomPairs': fputils.BuildAtomPairFP,
    'TopologicalTorsions': fputils.BuildTorsionsFP,
    'Avalon': fputils.BuildAvalonFP,
    'Gobbi2D': get_gobbi_fp,
    'MACCS': get_maccs_keys,
}

coeff_fn_dict = {
    'AllBit': DataStructs.AllBitSimilarity,
    'Asymmetric': DataStructs.AsymmetricSimilarity,
    'BraunBlanquet': DataStructs.BraunBlanquetSimilarity,
    'Cosine': DataStructs.CosineSimilarity,
    'Dice': DataStructs.DiceSimilarity,
    'Kulczynski': DataStructs.KulczynskiSimilarity,
    'McConnaughey': DataStructs.McConnaugheySimilarity,
    'OnBit': DataStructs.OnBitSimilarity,
    'RogotGoldberg': DataStructs.RogotGoldbergSimilarity,
    'Russel': DataStructs.RusselSimilarity,
    'Sokal': DataStructs.SokalSimilarity,
    'Tanimoto': DataStructs.TanimotoSimilarity,
    # 'Tversky': DataStructs.TverskySimilarity,
}

_FINGERPRINT_FN = fp_fn_dict['Morgan']  # default, app can override
_COEFFICIENT_FN = coeff_fn_dict['Tanimoto']  # default, app can override


def set_fingerprint_fn(fingerprint_name):
    global _FINGERPRINT_FN
    _FINGERPRINT_FN = fp_fn_dict[fingerprint_name]


def set_coefficient_fn(coeff_name):
    global _COEFFICIENT_FN
    _COEFFICIENT_FN = coeff_fn_dict[coeff_name]


def calculate_fingerprint(mol):
    return _FINGERPRINT_FN(mol)


def calculate_similarity(fp1, fp2):
    return _COEFFICIENT_FN(fp1, fp2)
