from datetime import datetime

from flask import Flask, g
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

    return app
