from time import time
from datetime import datetime
from flask import current_app
from flask_login import UserMixin

from . import db
from .. import paths
from ..admin import _get_md5


class User(UserMixin, db.Model):
    __tablename__ = 'users'
    id = db.Column(db.Integer, primary_key=True)
    social_id = db.Column(db.String(64), nullable=False, unique=True)
    display_name = db.Column(db.String(64), nullable=False)
    email = db.Column(db.String(64), nullable=True)
    last_seen = db.Column(db.DateTime, default=datetime.utcnow)
    in_team = db.Column(db.Boolean, default=False)
    alt_email_str = db.Column(db.String(255), nullable=True)
    is_admin = db.Column(db.Boolean, default=False)
    rebuilds = db.relationship('Rebuild', backref='user', lazy='dynamic')

    def __repr__(self):
        return '<User {}>'.format(self.email)

    @property
    def alt_emails(self):
        """Get list of non-primary emails associated with google ID."""
        if not self.alt_email_str:
            return []
        return self.alt_email_str.split(',')

    @property
    def known_emails(self):
        """Get set of all known emails for this user."""
        return {self.email}.union(set(self.alt_emails))

    def get_rebuilds_in_progress(self):
        return Rebuild.query.filter_by(user=self, complete=False).all()


class Rebuild(db.Model):
    __tablename__ = 'builds'
    id = db.Column(db.Integer, primary_key=True)
    start_time = db.Column(db.DateTime, index=True, default=datetime.utcnow)
    end_time = db.Column(db.DateTime, index=True, default=None)
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'))
    complete = db.Column(db.Boolean, default=False, index=True)
    status = db.Column(db.String(255))

    def __repr__(self):
        return '<Rebuild {}>'.format(self.id)

    def get_progress_message(self):
        # initial = f"Rebuild in progress."
        # return f"{initial} Started at {self.start_time}."
        return self.status

    def set_status_and_commit(self, message):
        self.status = message
        db.session.add(self)
        db.session.commit()

    def mark_complete_and_commit(self):
        self.complete = True
        self.end_time = datetime.utcnow()
        db.session.add(self)
        db.session.commit()

    @classmethod
    def get_rebuilds_in_progress(cls):
        return cls.query.filter_by(complete=False) \
            .order_by(Rebuild.start_time).all()

    @classmethod
    def get_most_recent_incomplete_rebuild(cls):
        return cls.query.filter_by(complete=False)\
            .order_by(Rebuild.start_time.desc()).first()

    @classmethod
    def get_most_recent_complete_rebuild(cls):
        return cls.query.filter_by(complete=True) \
            .order_by(Rebuild.end_time.desc()).first()


class ReferenceHash(db.Model):
    __tablename__ = 'reference_hash'
    is_gdrive = db.Column(db.Boolean, primary_key=True)
    hash = db.Column(db.String(255))

    @staticmethod
    def add_hash(use_drive=False, md5=None):
        version_str = 'gdrive' if use_drive else 'local'
        current_app.logger.debug(f'Setting {version_str} reference hash.')
        rh = ReferenceHash.query.get(use_drive)
        if rh is None:
            rh = ReferenceHash(is_gdrive=use_drive)
        rh.hash = md5
        db.session.add(rh)
        db.session.commit()

    @staticmethod
    def get_latest_hash_from_db(app):
        with app.app_context():
            use_drive = app.config['USE_DRIVE']
            rh = ReferenceHash.query.get(use_drive)
        md5 = rh.hash if rh is not None else None
        return md5

    # @staticmethod
    # def calculate_hash():
    #     md5 = _get_md5(paths.REFERENCE_PATH)
    #     return md5

    @staticmethod
    def update_and_get_hash():
        md5 = _get_md5(paths.REFERENCE_PATH)
        use_drive = current_app.config['USE_DRIVE']
        ReferenceHash.add_hash(use_drive=use_drive, md5=md5)
        return md5

    def __repr__(self):
        version = 'gdrive' if self.is_gdrive else 'local'
        return f'<ReferenceHash {version}>'
