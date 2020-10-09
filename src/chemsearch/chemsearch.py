"""
Run app with `flask run`
or
`gunicorn -b 0.0.0.0:5005 -w 2 --threads=4 --worker-class=gthread src.chemsearch.chemsearch:app`
"""
import os
import sys
import logging
import pathlib

# LOGGING
LOG_LEVEL = getattr(logging, os.environ.get('LOG_LEVEL', 'INFO').upper())
logging.basicConfig(format='%(levelname)s: %(message)s',  # %(asctime)-15s
                    level=LOG_LEVEL, stream=sys.stdout)
_logger = logging.getLogger(__name__)

# ENSURE DOTENV VARIABLES HAVE LOADED (for gunicorn)
from dotenv import load_dotenv
_env_path = str(pathlib.Path(__file__).parent.parent.parent
                .joinpath('config').joinpath('.env').absolute())
load_dotenv(_env_path, verbose=True)

from .app import create_app

app = create_app(os.getenv('FLASK_CONFIG') or 'default')


@app.shell_context_processor
def make_shell_context():
    from .app import db
    from .app.models import User, Rebuild
    return dict(db=db, User=User, Rebuild=Rebuild)


@app.cli.command()
def build():
    """Run deployment tasks."""
    from .app.rebuild import run_full_scan_and_rebuild
    with app.app_context():
        run_full_scan_and_rebuild(run_async=False)


@app.cli.command()
def link_data():
    from .app import link_data
    link_data(app)


@app.cli.command()
def deploy():
    """Run deployment tasks."""
    import shutil
    # TODO: create symlink in static to local_db dir
    use_env = '.env.deploy'
    shutil.copy(use_env, '.env')
    print(f"Copied from {use_env} to .env")


@app.cli.command()
def develop():
    """Set up development server."""
    import shutil
    # TODO: create symlink in static to local_db dir
    use_env = '.env.dev'
    shutil.copy(use_env, '.env')
    print(f"Copied from {use_env} to .env")
