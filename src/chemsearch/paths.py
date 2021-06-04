import os
import pathlib
from collections import defaultdict

from appdirs import AppDirs


def update_paths(use_drive=False):
    global ARCHIVE_DIR, PATH_SOURCES, SCAN_RESULTS_PATH, REFERENCE_PATH
    _version = 'gdrive' if use_drive else 'local'
    ARCHIVE_DIR = os.environ.get('LOCAL_DB_PATH')
    os.makedirs(ARCHIVE_DIR, exist_ok=True)
    PATH_SOURCES = _build_path_sources_dict(ARCHIVE_DIR)
    SCAN_RESULTS_PATH = PATH_SOURCES[_version]['scan']
    REFERENCE_PATH = PATH_SOURCES[_version]['master']


def _build_path_sources_dict(db_dir):
    sources = defaultdict(dict)
    for version in ['local', 'gdrive']:
        for source in ['scan', 'master']:
            sources[version][source] = os.path.join(db_dir,
                                                    f'{source}_{version}.tsv')
    return sources


_app_dirs = AppDirs('chemsearch')
DATA_ROOT = _app_dirs.user_data_dir
DEMO_DIR = pathlib.Path(_app_dirs.user_data_dir).joinpath('demo_db')
CONFIG_DIR = pathlib.Path(_app_dirs.user_config_dir).joinpath('config')
ENV_PATH = CONFIG_DIR.joinpath('.env')
SERVICE_ACCOUNT_CREDS = CONFIG_DIR.joinpath('creds.json')
SHORTCUTS_YAML = CONFIG_DIR.joinpath('custom_queries.yaml')
# DYNAMIC PATHS (set via update_paths)
ARCHIVE_DIR = None
PATH_SOURCES = None
SCAN_RESULTS_PATH = None
REFERENCE_PATH = None
