import os
import pathlib

ARCHIVE_DIR = os.environ.get('LOCAL_DB_PATH')
CONFIG_DIR = str(pathlib.Path(__file__).parent.parent.parent\
                 .joinpath('config').absolute())
ENV_PATH = os.path.join(CONFIG_DIR, '.env')
SCAN_RESULTS_PATH = os.path.join(ARCHIVE_DIR, 'scan_local.tsv')
REFERENCE_PATH = os.path.join(ARCHIVE_DIR, 'master_local.tsv')
SERVICE_ACCOUNT_CREDS = os.path.join(CONFIG_DIR, 'creds.json')


def update_paths(use_drive=False):
    global ARCHIVE_DIR, SCAN_RESULTS_PATH, REFERENCE_PATH
    ARCHIVE_DIR = os.environ.get('LOCAL_DB_PATH')
    if use_drive:
        SCAN_RESULTS_PATH = os.path.join(ARCHIVE_DIR, 'scan_gdrive.tsv')
        REFERENCE_PATH = os.path.join(ARCHIVE_DIR, 'master_gdrive.tsv')
    else:
        SCAN_RESULTS_PATH = os.path.join(ARCHIVE_DIR, 'scan_local.tsv')
        REFERENCE_PATH = os.path.join(ARCHIVE_DIR, 'master_local.tsv')
