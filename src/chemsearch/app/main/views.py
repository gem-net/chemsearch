import logging

from flask import render_template, flash, redirect, url_for, request, g, \
    current_app, abort
from flask_login import login_user, logout_user,\
    current_user

from . import main
from ..decorators import membership_required, admin_required
from .forms import admin_form_from_users, EmptyForm
from .. import db, META
from ..models import User
from .. import rebuild
from .. import filters
from ..oauth import OAuthSignIn
from ..paging import get_page_items_or_404, get_page_count
from ...db import get_substructure_matches, get_sim_matches, MolException, \
    LOCAL_MOLECULES, MOLECULE_DICT
from ... import drive


_logger = logging.getLogger(__name__)


@main.route('/', methods=['GET', 'POST'])
def index():
    if not current_app.config['USE_AUTH'] or not current_user.is_anonymous and current_user.in_team:
        pass_mols, filter_dict = filters.filter_mols(LOCAL_MOLECULES, request.args)
        filterable = filters.count_filterable(pass_mols)
        page_no = request.args.get('page', 1, type=int)
        molecules = get_page_items_or_404(pass_mols, page_no)
        n_pages = get_page_count(len(pass_mols))
        return render_template('index.html', molecules=molecules, n_pages=n_pages,
                               filters=filter_dict, filterable=filterable)
    else:
        return render_template('index.html')


@main.route('/molecule/<inchi_key>', methods=['GET', 'POST'])
@membership_required
def molecule(inchi_key):
    mol = MOLECULE_DICT.get(inchi_key)
    if mol is None:
        abort(404, 'Molecule not found.')
    return render_template('molecule.html', mol=mol)


@main.route('/custom/<inchi_key>', methods=['POST'])
@membership_required
def custom_info(inchi_key):
    mol = MOLECULE_DICT.get(inchi_key)
    if mol is None:
        abort(404, 'Molecule not found.')
    df, rec, content = drive.get_file_listing_and_custom_info(
        mol, files_resource=META.files_resource)
    custom_html = render_template('_custom.html', content=content, rec=rec)
    listing_html = render_template('_folder_files.html', df=df)
    return {'listing': listing_html, 'custom': custom_html}


@main.route('/search', methods=['GET', 'POST'])
@membership_required
def search():
    return render_template('search.html')


@main.route('/results', methods=['GET'])
@membership_required
def results():
    """
    Request args:
        smiles: str
        search_type: in ('similarity', 'substructure')
    """
    pass_mols, filter_dict = filters.filter_mols(LOCAL_MOLECULES, request.args)
    smiles = request.args.get('smiles')
    search_type = request.args.get('search_type')
    if smiles in (None, '') or search_type in (None, ''):
        flash("Bad inputs", "error")
        return redirect(url_for('.search'))
    page_no = request.args.get('page', 1, type=int)
    if search_type == 'substructure':
        try:
            molecules_all = get_substructure_matches(smiles, mols=pass_mols)
        except MolException as e:
            flash(str(e), "error")
            return redirect(url_for('.search'))
        n_pages = get_page_count(len(molecules_all))
        molecules = get_page_items_or_404(molecules_all, page_no) \
            if molecules_all else []
        sims = None  # no similarity scores for this search type
    elif search_type == 'similarity':
        try:
            sims_all, molecules_all = get_sim_matches(smiles, mols=pass_mols)
        except MolException as e:
            flash(str(e), "error")
            return redirect(url_for('.search'))
        n_pages = get_page_count(len(molecules_all))
        molecules = get_page_items_or_404(molecules_all, page_no)
        sims = get_page_items_or_404(sims_all, page_no)
    else:
        abort(404, "Unrecognized search type.")
    filterable = filters.count_filterable(molecules)
    return render_template('results.html', smiles=smiles,
                           molecules=molecules, sims=sims,
                           search_type=search_type, n_pages=n_pages,
                           filters=filter_dict, filterable=filterable)


@main.route('/admin', methods=['GET', 'POST'])
@admin_required
def admin():
    if current_app.config['USE_AUTH'] and not current_user.is_admin:
        abort(503, 'Unauthorized')
    if current_app.config['USE_AUTH']:
        other_users = User.query.filter(User.id != current_user.id).all()
        form = admin_form_from_users(current_user, other_users)
        if form.validate_on_submit():
            updated_users = form.update_admins(other_users)
            flash(f'{len(updated_users)} admins modified.')
            return redirect(url_for('main.admin'))
    else:
        form = None
    last_rebuild = rebuild.get_most_recent_complete_rebuild()
    in_progress_builds = rebuild.get_rebuilds_in_progress()
    empty_form = EmptyForm()
    return render_template('admin.html', user_form=form,
                           empty_form=empty_form,
                           last_rebuild=last_rebuild,
                           in_progress_builds=in_progress_builds,
                           )


@main.route('/clear-rebuilds', methods=['POST'])
@admin_required
def clear_rebuilds():
    in_progress = rebuild.get_rebuilds_in_progress()
    for r in in_progress:
        r.complete = None
        db.session.add(r)
    db.session.commit()
    flash("Incomplete rebuilds cleared.")
    return redirect(url_for('main.admin'))


@main.route('/logout')
def logout():
    logout_user()
    return redirect(url_for('.index'))


@main.route('/reload', methods=['GET', 'POST'])
@admin_required
def reload():
    from ..rebuild import run_full_scan_and_rebuild
    run_full_scan_and_rebuild(current_user)
    return redirect(url_for('main.admin'))


@main.route('/authorize/<provider>')
def oauth_authorize(provider):
    if not current_user.is_anonymous:
        return redirect(url_for('.index'))
    oauth_obj = OAuthSignIn.get_provider(provider)
    return oauth_obj.authorize()


@main.route('/callback/<provider>')
def oauth_callback(provider):
    if not current_user.is_anonymous:
        return redirect(url_for('.index'))
    oauth_obj = OAuthSignIn.get_provider(provider)
    social_id, username, email, alt_email_str = oauth_obj.callback()
    if social_id is None:
        flash('Authentication failed.', 'error')
        return redirect(url_for('.index'))
    user = User.query.filter_by(social_id=social_id).first()
    modified = False
    # CREATE NEW USER IF NOT IN DB
    if not user:
        if social_id in g.members_dict:
            email = g.members_dict[social_id]
        user = User(social_id=social_id, display_name=username, email=email,
                    alt_email_str=alt_email_str)
        _logger.info(f"Creating new user: {email}.")
        modified = True
    # UPDATE EMAIL DATA IF NECESSARY
    if user.email != email or user.alt_email_str != alt_email_str:
        user.email = email
        user.alt_email_str = alt_email_str
        _logger.info(f"Modifying email for user: {email}.")
        modified = True
    if modified:
        db.session.add(user)
        db.session.commit()
    login_user(user, True)
    _logger.info(f"Logging in user: {email}.")
    return redirect(url_for('.index'))
