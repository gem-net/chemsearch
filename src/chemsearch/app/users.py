import os
import logging
from collections import namedtuple

from flask import g, current_app
from google.oauth2 import service_account
from googleapiclient.discovery import build

from . import login_manager
from .models import User

logger = logging.getLogger(__name__)
MEMBERS_DICT = {}  # updated at startup with update_members_dict_from_config
DIR_SERVICE_HANDLE = None  # set at app startup


@login_manager.user_loader
def load_user(user_id):
    if not current_app.config['USE_AUTH']:
        return None
    return User.query.get(int(user_id))


def _set_service_handle_using_config(app):
    """Get dictionary of {service_name: service_handle}, else return None."""
    global DIR_SERVICE_HANDLE
    if not app.config['USE_DRIVE'] and not app.config['USE_AUTH']:
        app.logger.info("Skipping DIR service creation (Drive/Auth only).")
        return
    service_account_file = app.config['SERVICE_ACCOUNT_FILE']
    scopes = [
        'https://www.googleapis.com/auth/admin.directory.user.readonly',
        'https://www.googleapis.com/auth/admin.directory.group.member.readonly',
        ]
    credentials = service_account.Credentials.from_service_account_file(
        service_account_file, scopes=scopes)
    delegated_credentials = credentials.with_subject(
        app.config['CREDENTIALS_AS_USER'])
    DIR_SERVICE_HANDLE = build('admin', 'directory_v1',
                               credentials=delegated_credentials,
                               cache_discovery=False)


def _get_members_dict(app):
    """Get dictionary of {google_id: email_address}.

    Return:
         members_dict (dict): google ID: email dictionary.
    """
    service = DIR_SERVICE_HANDLE

    # GOOGLE GROUP MEMBERS
    group_key = app.config.get('GROUP_KEY', None)
    logging.info("Looking up group members list.")
    res = service.members().list(groupKey=group_key).execute()
    members_dict = {i['id']: i['email'] for i in res['members'] if 'email' in i}

    # DOMAIN-ONLY MEMBERS
    domain_users = _get_domain_users()
    domain_dict = {u.id: u.email for u in domain_users}
    members_dict.update(domain_dict)

    return members_dict


def _get_domain_users():
    """Get list of domain users (who might not be in specified Google Group).

    Warning: will only fetch up to 100 users.
    """
    users = []
    service = DIR_SERVICE_HANDLE
    logger.info("Looking up domain-specific users.")
    res = service.users().list(customer='my_customer').execute()
    UserInfo = namedtuple('UserInfo', ['id', 'email', 'full_name'])
    for user in res['users']:
        user_id = user['id']
        full_name = user['name']['fullName']
        email = user['primaryEmail']
        users.append(UserInfo(user_id, email, full_name))
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


def update_members_dict_from_config(app):
    """Update members dictionary."""
    if not app.config['USE_AUTH']:
        return
    if DIR_SERVICE_HANDLE is None:
        _set_service_handle_using_config(app)
    new_dict = _get_members_dict(app)
    MEMBERS_DICT.clear()
    MEMBERS_DICT.update(new_dict)
