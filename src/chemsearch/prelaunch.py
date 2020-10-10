"""
For creating initial db in Docker container.
"""
import os

from .app import db, Flask
from .app.config import config
from .paths import ENV_PATH

# ENSURE DOTENV VARIABLES HAVE LOADED (for gunicorn)
from dotenv import load_dotenv
load_dotenv(ENV_PATH, verbose=True)

config_name = os.getenv('FLASK_CONFIG') or 'default'
app = Flask(__name__)
app.config.from_object(config[config_name])
db.init_app(app)

from .app.models import User, Rebuild

with app.app_context():
    db.create_all()

print("DB create_all complete.")
