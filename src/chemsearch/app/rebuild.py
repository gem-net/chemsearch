from . import db
from .models import Rebuild


def get_rebuilds_in_progress():
    return Rebuild.query.filter_by(complete=False).all()


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
