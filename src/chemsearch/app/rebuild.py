import os
import logging
from threading import Thread

from flask import current_app

from . import db
from .models import Rebuild

from .. import drive
from ..admin import assemble_archive_metadata
from ..db import reload_molecules


def run_full_scan_and_rebuild(build_id):
    """Rescan and rebuild local archive. Requires app context."""
    app = current_app._get_current_object()
    thr = Thread(target=run_full_scan_and_rebuild_async, args=[app, build_id])
    thr.start()
    return thr


def run_full_scan_and_rebuild_async(app, build_id: str):
    from .. import _logger
    with app.app_context():
        archive_dir = current_app.config['LOCAL_DB_PATH']
        log_path = os.path.join(archive_dir, f'rebuild_{build_id}.log')
        fh = logging.FileHandler(log_path)
        fh.setLevel(logging.DEBUG)
        _logger.addHandler(fh)

        build = Rebuild.query.get(build_id)  # type: Rebuild
        build.set_status_and_commit("Identifying categories and MOL files in Drive")
        meta = drive.Meta().build()

        build.set_status_and_commit("Updating local archive.")
        drive.create_local_archive(meta.molfiles, local_root=archive_dir,
                                   files_resource=meta.files_resource)
        build.set_status_and_commit("Generating images and metadata.")
        df = assemble_archive_metadata(archive_dir)
        build.set_status_and_commit(f"Completed rebuild contains {len(df)} molecules.")
        build.mark_complete_and_commit()

        reload_molecules()
        _logger.removeHandler(fh)


def get_rebuilds_in_progress():
    return Rebuild.query.filter_by(complete=False)\
        .order_by(Rebuild.start_time).all()


def get_most_recent_complete_rebuild():
    return Rebuild.query.filter_by(complete=True)\
        .order_by(Rebuild.end_time.desc()).first()


def mark_rebuilds_as_failed(rebuild_list, commit=True):
    for rebuild in rebuild_list:
        rebuild.complete = None
        if commit:
            db.session.add(rebuild)
    if commit:
        db.session.commit()
