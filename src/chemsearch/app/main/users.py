import os
import logging
from datetime import datetime
from collections import namedtuple

from flask import g
from flask_login import UserMixin
from google.oauth2 import service_account
from googleapiclient.discovery import build

from .. import db, login_manager


logger = logging.getLogger(__name__)
MEMBERS_DICT = {}  # updated at end of module
# DIR_SERVICE_HANDLE = None  # saved at end of module


class User(UserMixin, db.Model):
    __tablename__ = 'users'
    id = db.Column(db.Integer, primary_key=True)
    social_id = db.Column(db.String(64), nullable=False, unique=True)
    display_name = db.Column(db.String(64), nullable=False)
    email = db.Column(db.String(64), nullable=True)
    last_seen = db.Column(db.DateTime, default=datetime.utcnow)
    in_cgem = db.Column(db.Boolean, default=False)
    alt_email_str = db.Column(db.String(255), nullable=True)

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


@login_manager.user_loader
def load_user(user_id):
    return User.query.get(int(user_id))


def _get_dir_service_handle():
    """Get dictionary of {service_name: service_handle}."""
    service_account_file = os.getenv('SERVICE_ACCOUNT_FILE')
    scopes = [
        'https://www.googleapis.com/auth/admin.directory.user.readonly',
        'https://www.googleapis.com/auth/admin.directory.group.member.readonly',
        ]

    credentials = service_account.Credentials.from_service_account_file(
        service_account_file, scopes=scopes)
    delegated_credentials = credentials.with_subject(
        os.environ.get('CREDENTIALS_AS_USER'))
    dir_service = build('admin', 'directory_v1', credentials=delegated_credentials,
                        cache_discovery=False)

    return dir_service


def get_members_dict():
    """Get dictionary of {google_id: email_address}.

    Return:
         members_dict (dict): google ID: email dictionary.
    """
    service = DIR_SERVICE_HANDLE

    # GOOGLE GROUP MEMBERS
    group_key = os.environ.get('GROUP_KEY', None)
    logging.info("Looking up group members list.")
    res = service.members().list(groupKey=group_key).execute()
    members_dict = {i['id']: i['email'] for i in res['members'] if 'email' in i}

    # DOMAIN-ONLY MEMBERS
    domain_users = get_domain_users()
    domain_dict = {u.id: u.email for u in domain_users}
    members_dict.update(domain_dict)

    return members_dict


def get_domain_users():
    """Get list of domain users (who might not be in specified Google Group).

    Warning: will only fetch up to 100 users.
    """
    users = []
    service = DIR_SERVICE_HANDLE
    logger.info("Looking up domain-specific users.")
    res = service.users().list(customer='my_customer').execute()
    User = namedtuple('User', ['id', 'email', 'full_name'])
    for user in res['users']:
        user_id = user['id']
        full_name = user['name']['fullName']
        email = user['primaryEmail']
        users.append(User(user_id, email, full_name))
    return users


def find_user_by_email(find_email):
    user = None
    for u in User.query.all():
        if find_email in u.known_emails:
            user = u
            break
    return user


def update_g():
    g.members_dict = MEMBERS_DICT


def update_members_dict():
    """Update members dictionary"""
    global MEMBERS_DICT
    MEMBERS_DICT.clear()
    new_dict = get_members_dict()
    MEMBERS_DICT.update(new_dict)


DIR_SERVICE_HANDLE = _get_dir_service_handle()
update_members_dict()
