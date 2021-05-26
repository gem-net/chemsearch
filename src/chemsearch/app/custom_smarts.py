import hashlib
import logging
from collections import namedtuple

from rdkit import Chem
from flask import Flask

from . import db
from .models import CustomSpec, CustomMatch
from ..db import LOCAL_MOLECULES


_logger = logging.getLogger(__name__)

Spec = namedtuple('Spec', ['query_name', 'smarts_str', 'digest'])

SPEC_LIST = []  # type: list[Spec] # initialized from app config with init_app


def lookup_custom_matches(app: Flask, spec: Spec):
    with app.app_context():
        match_mols = []
        for lm in LOCAL_MOLECULES:
            cm = CustomMatch.query.get((spec.digest, lm.inchi_key))
            if cm.sub_match:
                match_mols.append(lm)
    return match_mols


def update_custom_spec_db(app: Flask):
    with app.app_context():
        _remove_discarded_db_specs(SPEC_LIST)
        _add_db_specs(SPEC_LIST)
        _fill_custom_spec_matches()


def _load_custom_spec_list(app):
    """Gather custom query Spec tuples."""
    custom_dict = app.config['CUSTOM_QUERIES']
    spec_list = []
    for query_name, smarts_str in custom_dict.items():
        # Report failed mol construction
        query_mol = Chem.MolFromSmarts(smarts_str)
        if query_mol is None:
            _logger.info(f"Failed to build Mol for custom query: {query_name}. Skipping.")
            continue
        digest = _get_digest(smarts_str)
        spec = Spec(query_name=query_name, smarts_str=smarts_str, digest=digest)
        spec_list.append(spec)
    return spec_list


def _remove_discarded_db_specs(spec_list):
    """Requires app context."""
    digests = {i.digest for i in spec_list}
    specs_db = CustomSpec.query.all()
    is_changed = False
    for rec in specs_db:
        if rec.smarts_hash not in digests:
            is_changed = True
            _logger.info(f"Removing old query with unrecognised SMARTS: {rec.query_name}.")
            db.session.delete(rec)
    if is_changed:
        db.session.commit()
    else:
        _logger.info("No discarded custom queries in db.")


def _add_db_specs(spec_list):
    """Requires app context."""
    added_smarts = False
    for spec in spec_list:
        cs = CustomSpec.query.get(spec.digest)
        # create smarts row if not found
        if cs is None:
            _logger.info(f"Adding custom SMARTS to db: {spec.query_name}.")
            cs = CustomSpec(smarts_hash=spec.digest, query_name=spec.query_name,
                            smarts_str=spec.smarts_str)
            db.session.add(cs)
            added_smarts = True
        # modify name if mismatch
        elif spec.query_name != cs.query_name:
            _logger.info(f"Renaming custom SMARTS: {cs.query_name} -> {spec.query_name}.")
            cs.query_name = spec.query_name
            db.session.add(cs)
            modified = True
    db.session.commit()
    if not added_smarts:
        _logger.info("No new custom queries to add.")


def _fill_custom_spec_matches():
    cs_list = CustomSpec.query.all()
    n_added = 0
    for lm in LOCAL_MOLECULES:
        for cs in cs_list:
            # add match info if not already present
            if CustomMatch.query.get((cs.smarts_hash, lm.inchi_key)) is None:
                sub_match = lm.has_substructure(cs.smarts_str)
                cm = CustomMatch(smarts_hash=cs.smarts_hash, inchi_key=lm.inchi_key, sub_match=sub_match)
                db.session.add(cm)
                n_added += 1
    if n_added:
        _logger.info(f"Adding match info for {n_added} MOL-SMARTS pairs.")
        db.session.commit()
    else:
        _logger.info("No new match info to add to db.")


def _get_digest(val):
    return hashlib.md5(val.encode('utf8')).hexdigest()


def init_app(app):
    global SPEC_LIST
    SPEC_LIST.clear()
    SPEC_LIST.extend(_load_custom_spec_list(app))
