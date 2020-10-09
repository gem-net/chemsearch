from typing import List

from flask_wtf import FlaskForm
from wtforms import BooleanField, SubmitField

from .. import db
from ..users import User


class AdminForm(FlaskForm):
    """Customized via admin_form_from_users."""
    def update_admins(self, other_users):
        updated_users = []
        for user in other_users:
            user_key = f"user_{user.id}"
            if user_key in self.data:
                marked_admin = self.data[user_key]
                if marked_admin != user.is_admin:
                    user.is_admin = marked_admin
                    db.session.add(user)
                    updated_users.append(user)
        if updated_users:
            db.session.commit()
        return updated_users


def admin_form_from_users(current_user: User, other_users: List[User]):
    class CustomAdminForm(AdminForm):
        pass
    if not current_user.is_anonymous:
        field = BooleanField(label=current_user.email, default="checked")
        setattr(CustomAdminForm, 'current_user', field)
    for user in other_users:
        field_name = f"user_{user.id}"
        is_admin = user.is_admin
        default = "checked" if is_admin else None
        field = BooleanField(label=user.email, default=default)
        setattr(CustomAdminForm, field_name, field)
    CustomAdminForm.submit = SubmitField(label='Update')

    form = CustomAdminForm()
    return form


class EmptyForm(FlaskForm):
    submit = SubmitField('Submit')
