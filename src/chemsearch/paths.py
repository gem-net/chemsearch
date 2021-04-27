import os
import pathlib
from collections import defaultdict


def update_paths(use_drive=False):
    global ARCHIVE_DIR, PATH_SOURCES, SCAN_RESULTS_PATH, REFERENCE_PATH
    _version = 'gdrive' if use_drive else 'local'
    ARCHIVE_DIR = os.environ.get('LOCAL_DB_PATH')
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


ROOT_DIR = pathlib.Path(__file__).parent.parent.parent
DEMO_DIR = ROOT_DIR.joinpath('demo_db')
CONFIG_DIR = ROOT_DIR.joinpath('config').absolute()
ENV_PATH = CONFIG_DIR.joinpath('.env')
SERVICE_ACCOUNT_CREDS = CONFIG_DIR.joinpath('creds.json')
ARCHIVE_DIR = None
PATH_SOURCES = None
SCAN_RESULTS_PATH = None
REFERENCE_PATH = None
