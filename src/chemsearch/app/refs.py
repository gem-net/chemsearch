import os
import yaml

from ..paths import CONFIG_DIR

TEMPLATES = dict()


def update_external_db_templates():
    yaml_path = os.path.join(CONFIG_DIR, 'external_dbs.yaml')
    if os.path.exists(yaml_path):
        with open(yaml_path, 'r') as infile:
            templates = yaml.load(infile, Loader=yaml.SafeLoader)
        TEMPLATES.clear()
        TEMPLATES.update(templates)


update_external_db_templates()
