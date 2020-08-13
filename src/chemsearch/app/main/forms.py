from typing import List

from flask_wtf import FlaskForm
from wtforms import BooleanField, SubmitField

from ..users import User


def admin_form_from_users(current_user: User, other_users: List[User]):
    class AdminForm(FlaskForm):
        pass
    if not current_user.is_anonymous:
        field = BooleanField(label=current_user.email, default="checked")
        setattr(AdminForm, 'current_user', field)
    for user in other_users:
        field_name = f"user_{user.id}"
        is_admin = user.is_admin
        default = "checked" if is_admin else None
        field = BooleanField(label=user.email, default=default)
        setattr(AdminForm, field_name, field)
        AdminForm.submit = SubmitField(label='Update')

    form = AdminForm()
    return form
