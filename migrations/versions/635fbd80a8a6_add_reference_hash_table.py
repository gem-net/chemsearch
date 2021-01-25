"""Add reference_hash table.

Revision ID: 635fbd80a8a6
Revises: 7d9f1fe19b3b
Create Date: 2021-01-20 13:26:30.454630

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '635fbd80a8a6'
down_revision = '7d9f1fe19b3b'
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.create_table('reference_hash',
    sa.Column('is_gdrive', sa.Boolean(), nullable=False),
    sa.Column('hash', sa.String(length=255), nullable=True),
    sa.PrimaryKeyConstraint('is_gdrive')
    )
    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_table('reference_hash')
    # ### end Alembic commands ###
