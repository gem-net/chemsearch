import os

ARCHIVE_DIR = os.environ.get('LOCAL_DB_PATH')
SCAN_RESULTS_PATH = os.path.join(ARCHIVE_DIR, 'scan_local.tsv')
REFERENCE_PATH = os.path.join(ARCHIVE_DIR, 'master_local.tsv')


def update_paths(use_drive=False):
    global ARCHIVE_DIR, SCAN_RESULTS_PATH, REFERENCE_PATH
    ARCHIVE_DIR = os.environ.get('LOCAL_DB_PATH')
    if use_drive:
        SCAN_RESULTS_PATH = os.path.join(ARCHIVE_DIR, 'scan_gdrive.tsv')
        REFERENCE_PATH = os.path.join(ARCHIVE_DIR, 'master_gdrive.tsv')
    else:
        SCAN_RESULTS_PATH = os.path.join(ARCHIVE_DIR, 'scan_local.tsv')
        REFERENCE_PATH = os.path.join(ARCHIVE_DIR, 'master_local.tsv')
