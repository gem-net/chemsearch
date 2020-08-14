from time import time
from datetime import datetime
from flask_login import UserMixin

from . import db


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

    def __repr__(self):
        return '<Rebuild {}>'.format(self.id)

    def get_progress_message(self):
        initial = f"Rebuild in progress."
        return f"{initial} Started at {self.start_time}."


# class Status(db.Model):
#     __tablename__ = 'status'
#     last_mol_modtime = db.Column(db.DateTime, index=True,
#                                  default=datetime.fromtimestamp(0))
#     last_refresh_time = db.Column(db.DateTime, index=True,
#                                   default=datetime.fromtimestamp(0))

