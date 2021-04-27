"""
Run app with `flask run`
or
`gunicorn -b 0.0.0.0:5005 -w 2 --threads=4 --worker-class=gthread src.chemsearch.chemsearch:app`
"""
import os
import re
import click
import logging
import shutil
import pathlib
from typing import Union

import dotenv
from flask.cli import FlaskGroup

from . import logger
from .paths import ENV_PATH, CONFIG_DIR

# ENSURE DOTENV VARIABLES HAVE LOADED (for gunicorn)
DOTENV_VALS = dotenv.dotenv_values(ENV_PATH)
dotenv.load_dotenv(ENV_PATH, verbose=False)  # Allow missing .env
if 'FLASK_APP' not in os.environ:
    os.environ['FLASK_APP'] = 'chemsearch.chemsearch'

from .app import create_app, link_data, init_data

app = create_app(os.getenv('FLASK_ENV') or 'default')


# REFERENCE
# @click.option('-y', '--yaml', 'query_yaml', required=True,
#               type=click.Path(exists=True),
#               help='Path to YAML file with title and abstract fields.')


@click.group()
@click.option('--verbose/--quiet', default=False)
@click.pass_context
def cli(ctx, verbose):
    if verbose:
        logger.setLevel("DEBUG")
    elif ctx.invoked_subcommand != 'flask':
        logger.handlers[0].setFormatter(logging.Formatter('%(message)s'))


# @cli.command()
# def run(**kwargs):
#     """Run the server."""
#     app.run(port=os.environ.get('FLASK_RUN_PORT', 5000))


def create_app_cli():
    env_name = os.environ.get('FLASK_ENV', 'production')
    return create_app(env_name)


@cli.group(cls=FlaskGroup, create_app=create_app_cli)
def flask():
    """Run Flask server commands."""
    # TODO: ensure config env is set up.


@cli.group()
def setup():
    """Set up chemsearch, configuring ENV variables."""


@setup.command('config')
@click.option('-d/-l', '--use-drive/--local-only', default=None,
              show_default=True, help="Use DRIVE mode, linking Shared Drive.")
@click.option('-a/-n', '--use-auth/--no-auth', default=None, show_default=True,
              help="Use AUTH mode, with Google authentication.")
def configure(use_drive, use_auth):
    """Create configuration .env file via command prompt."""
    env_dict = dict()
    prev_vals = DOTENV_VALS.copy()

    def _set_env_via_prompt(var_name, prompt_str, default=None, **prompt_kw):
        if var_name in prev_vals:
            default = prev_vals[var_name]
        env_dict[var_name] = click.prompt(prompt_str, default=default, **prompt_kw)

    if use_drive is None:
        use_drive = click.confirm(
            "Use DRIVE mode (with Shared Drive storage)?", default=False)
    if use_auth is None:
        use_auth = click.confirm(
            "Use AUTH mode (with Google Group authentication)?", default=False)

    bold_style = dict(bold=True, fg='red')
    for test, mode_str, mode_desc in [
                (use_drive and use_auth, 'DRIVE+AUTH', 'DRIVE and AUTH mode'),
                (use_drive, 'DRIVE', 'DRIVE mode without authentication'),
                (use_auth, 'AUTH', 'LOCAL mode with authentication'),
                (True, 'BASIC', 'LOCAL mode without authentication')]:
        if test:
            click.secho(f"Building configuration for {mode_desc}.", **bold_style)
            break

    _set_env_via_prompt('LOCAL_DB_PATH', 'Local archive folder',
                        default=os.getcwd(),
                        type=click.Path(exists=True, file_okay=False,
                                        dir_okay=True, writable=True,
                                        resolve_path=True))

    # FLASK_RUN_PORT not used unless flask run from .env directory
    # _set_env_via_prompt('FLASK_RUN_PORT', 'Port', default=5000)
    _set_env_via_prompt('APP_TITLE', 'App Title', default='Chemsearch')
    _set_env_via_prompt('FLASK_ENV', 'Environment', default='production',
                        type=click.Choice(['development', 'production']))
    if use_drive or use_auth:
        _set_env_via_prompt('CREDENTIALS_AS_USER',
                            f'Authorized email (for {mode_str} credentials)',
                            type=EMAIL_TYPE)
    if use_drive:
        _set_env_via_prompt('SHARED_DRIVE_ID',
                            'Shared Drive ID (e.g. 0ABC1234DEF5)')
    if use_auth:
        _set_env_via_prompt('GOOGLE_CLIENT_ID', 'Google OAuth ClientID')
        _set_env_via_prompt('GOOGLE_SECRET', 'Google OAuth secret')
        _set_env_via_prompt('GROUP_KEY', 'Members Google Group ID')

    if ENV_PATH.exists():
        msg = "OK to overwrite previous .env file with these values?"
    else:
        msg = "OK to create an .env file with these values?"
    is_happy = click.confirm(msg, default=True)
    if not is_happy:
        click.echo("Aborting .env file creation.")
    backup_path = store_new_env(env_dict)
    if backup_path:
        click.echo(f"Previous config file saved to {backup_path}")

    build_ok = click.confirm("Build metadata?", default=True)
    if build_ok:
        init_data(app)
    else:
        click.echo("You can extract metadata later with 'rebuild' on Admin page.")


@setup.command()
def revert():
    """Revert to previous .env file."""
    prev_path = _get_previous_env_path()
    if prev_path.exists():
        overwrite = click.confirm("Are you sure you want to restore previous env file?")
        if overwrite:
            prev_path.rename(ENV_PATH)
            click.echo("Configuration reset to previous values.")
        else:
            click.echo("Aborted.")
    else:
        click.echo("No previous configuration file found.")


@setup.command()
def show():
    """Print configuration path and contents."""
    if ENV_PATH.exists():
        click.echo(f"ENV PATH: {ENV_PATH}\n\nContents:")
        with open(ENV_PATH, 'r') as env:
            for line in env:
                click.secho(line.strip(), fg='green')
    else:
        click.echo("No configuration file found.")


@setup.command()
def edit():
    """Open configuration .env file in an editor."""
    if not ENV_PATH.exists():
        ENV_PATH.touch()
    click.echo(f"Opening {ENV_PATH}. Edit, then save when you're done.")
    click.launch(str(ENV_PATH))


@setup.command()
@click.argument('path', type=click.Path(exists=True), required=True)
def creds(path):
    """Copy Google JSON credentials to config folder."""
    if not path.lower().endswith('.json'):
        click.echo("Please provide a path that ends with '.json'.")
        return
    shutil.copy(path, CONFIG_DIR.joinpath('creds.json'))
    click.echo("Credentials file copied successfully.")


@setup.command('import')
@click.argument('path', type=click.Path(exists=True), required=True)
def load(path):
    """Load variables from specified .env path."""
    click.echo(f"{path=}")


def store_new_env(env_dict) -> Union[None, os.PathLike]:
    backup_path = None
    if ENV_PATH.exists():
        backup_path = _get_previous_env_path()
        ENV_PATH.rename(backup_path)
    ENV_PATH.touch()
    for key in env_dict:
        dotenv.set_key(ENV_PATH, key, env_dict[key])
    logger.info("Created .env file")
    return backup_path


def _get_previous_env_path():
    backup_name = ENV_PATH.name + '.prev'
    return ENV_PATH.parent.joinpath(backup_name)


class EmailType(click.ParamType):
    name = "email"

    def convert(self, value, param, ctx):
        if re.fullmatch('[^@]+@[^@]+\.[^@]+', value):
            return value
        else:
            self.fail(f"{value!r} is not a valid email address.", param, ctx)


EMAIL_TYPE = EmailType()


@cli.command()
def build():
    """Run deployment tasks."""
    from .app.rebuild import run_full_scan_and_rebuild
    with app.app_context():
        run_full_scan_and_rebuild(run_async=False)


@cli.command()
def link_data():
    """Create data symlink in static directory."""
    link_data(app)


@cli.command()
def deploy():
    """Run deployment tasks."""
    import shutil
    # TODO: create symlink in static to local_db dir
    use_env = '.env.deploy'
    # shutil.copy(use_env, '.env')
    print(f"Mock: Copied from {use_env} to .env")


@cli.command()
def develop():
    """Set up development server."""
    import shutil
    # TODO: create symlink in static to local_db dir
    use_env = '.env.dev'
    shutil.copy(use_env, '.env')
    print(f"Copied from {use_env} to .env")


@app.shell_context_processor
def make_shell_context():
    from .app import db
    from .app.models import User, Rebuild
    return dict(db=db, User=User, Rebuild=Rebuild)
