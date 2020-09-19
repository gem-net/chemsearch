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


bootstrap = Bootstrap()
db = SQLAlchemy()
login_manager = LoginManager()
migrate = Migrate()
moment = Moment()
# login_manager.login_view = 'main.index'


def create_app(config_name):
    app = Flask(__name__)
    app.config.from_object(config[config_name])
    config[config_name].init_app(app)

    bootstrap.init_app(app)
    db.init_app(app)
    login_manager.init_app(app)
    migrate.init_app(app, db)
    moment.init_app(app)

    from ..admin import update_paths
    update_paths(use_drive=app.config['USE_DRIVE'])

    # session.init_app(app)

    from .main import main as main_blueprint
    app.register_blueprint(main_blueprint)
    # update_members_dict(app)

    with app.app_context():
        db.create_all()

    @app.before_request
    def before_request():
        from .users import update_g
        update_g()
        if app.config['USE_AUTH']:
            if current_user.is_authenticated:
                current_user.last_seen = datetime.utcnow()
                # first-request data reload might be necessary for auto-restart development
                current_user.in_team = current_user.social_id in g.members_dict
                db.session.commit()

    @app.before_first_request
    def before_first_request():
        link_data(app)

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
    if os.path.exists(data_path):
        current_target = os.readlink(data_path)
        if current_target != local_db_path:
            app.logger.info(f"Updating static/data target.")
            os.remove(data_path)
            symlink_needed = True
        else:
            app.logger.info(f"Current static/data target is correct.")
    else:
        symlink_needed = True
        app.logger.info("Creating static/data link.")
    if symlink_needed:
        data_path.symlink_to(local_db_path, target_is_directory=True)
        app.logger.info(f"New data directory is {local_db_path}.")
