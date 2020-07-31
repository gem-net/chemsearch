import logging

from flask import render_template, flash, redirect, url_for, request

from . import main
from ...db import get_substructure_matches, get_sim_matches, MolException

_logger = logging.getLogger(__name__)


@main.route('/', methods=['GET', 'POST'])
def index():
    from chemsearch.db import get_molecules
    molecules = list(get_molecules())

    return render_template('index.html', molecules=molecules)


@main.route('/search', methods=['GET', 'POST'])
def search():
    return render_template('search.html')


@main.route('/results', methods=['GET'])
def results():
    """
    Request args:
        smiles: str
        search_type: in ('similarity', 'substructure')
    """
    smiles = request.args.get('smiles')
    search_type = request.args.get('search_type')
    if smiles in (None, '') or search_type in (None, ''):
        flash("Bad inputs", "error")
        return redirect(url_for('.search'))
    if search_type == 'substructure':
        try:
            molecules = get_substructure_matches(smiles)
        except MolException as e:
            flash(str(e), "error")
            return redirect(url_for('.search'))
        sims = None  # no similarity scores for this search type
    else:
        try:
            sims, molecules = get_sim_matches(smiles)
        except MolException as e:
            flash(str(e), "error")
            return redirect(url_for('.search'))
    return render_template('results.html', smiles=smiles,
                           molecules=molecules,
                           sims=sims,
                           search_type=search_type)
