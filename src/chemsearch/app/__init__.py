import os
import pathlib
from datetime import datetime

from flask import Flask, g, request
from flask_bootstrap import Bootstrap
from flask_sqlalchemy import SQLAlchemy
from flask_login import current_user, LoginManager
from flask_migrate import Migrate
from flask_moment import Moment

# from flask_session import Session

from .config import config
from ..drive import Meta


bootstrap = Bootstrap()
db = SQLAlchemy(session_options={'expire_on_commit': False})
login_manager = LoginManager()
migrate = Migrate()
moment = Moment()
# login_manager.login_view = 'main.index'

META = Meta()


def create_app(config_name):
    app = Flask(__name__)
    app.config.from_object(config[config_name])
    config[config_name].init_app(app)

    bootstrap.init_app(app)
    db.init_app(app)
    login_manager.init_app(app)
    migrate.init_app(app, db)
    moment.init_app(app)

    from ..paths import update_paths
    update_paths(use_drive=app.config['USE_DRIVE'])

    from .main import main as main_blueprint
    app.register_blueprint(main_blueprint)

    # SET CONFIG_DEPENDENT MODULE ATTRIBUTES
    from .filters import set_filters_using_config
    set_filters_using_config(app)

    from . import rebuild
    from .models import ReferenceHash
    from ..db import reload_molecules

    @app.before_first_request
    def before_first_request():
        from .users import update_members_dict_from_config
        update_members_dict_from_config(app)
        link_data(app)
        rebuild.CURRENT_REF_HASH = ReferenceHash.get_latest_hash_from_db(app)
        reload_molecules()

    @app.before_request
    def before_request():
        from .users import MEMBERS_DICT
        g.members_dict = MEMBERS_DICT
        from .refs import TEMPLATES
        g.db_templates = TEMPLATES
        if app.config['USE_AUTH']:
            if current_user.is_authenticated:
                current_user.last_seen = datetime.utcnow()
                # first-request data reload might be necessary for auto-restart development
                current_user.in_team = current_user.social_id in g.members_dict
                db.session.commit()
        latest_hash = ReferenceHash.get_latest_hash_from_db(app)
        if rebuild.CURRENT_REF_HASH != latest_hash:
            app.logger.info("Identified updated reference hash "
                            f"({rebuild.CURRENT_REF_HASH} > {latest_hash}). "
                            "Reloading molecules.")
            from ..db import reload_molecules
            reload_molecules()
            rebuild.CURRENT_REF_HASH = latest_hash

    @app.context_processor
    def utility_processor():
        from .filters import update_args

        def get_updated_args(attr, val):
            return update_args(request.args, attr, val)

        return dict(get_updated_args=get_updated_args)

    return app


def link_data(app):
    local_db_path = os.path.abspath(app.config['LOCAL_DB_PATH'])
    static_path = app.static_folder
    data_path = pathlib.Path(static_path).joinpath('data')
    symlink_needed = False
    if data_path.is_symlink():
        current_target = os.readlink(data_path)
        if current_target != local_db_path:
            app.logger.info(f"Updating static/data target: {current_target} -> {local_db_path}.")
            try:
                os.remove(data_path)
            except FileNotFoundError:
                pass  # handle multiprocess race condition
            symlink_needed = True
        else:
            app.logger.debug(f"Current static/data target is correct.")
    else:
        app.logger.info("Creating static/data link.")
        symlink_needed = True
    if symlink_needed:
        try:
            data_path.symlink_to(local_db_path, target_is_directory=True)
            app.logger.info(f"New data directory is {local_db_path}.")
        except FileExistsError:  # handle multiprocess race condition
            app.logger.debug(f"Skipping overwrite of existing data dir.")
