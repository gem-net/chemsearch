import yaml
import pathlib

TEMPLATES = dict()


def update_external_db_templates():
    package_root = pathlib.Path(__file__).parent.parent.parent.parent
    yaml_path = package_root.joinpath('external_dbs.yaml')
    with open(yaml_path, 'r') as infile:
        templates = yaml.load(infile, Loader=yaml.SafeLoader)
    TEMPLATES.clear()
    TEMPLATES.update(templates)


update_external_db_templates()
