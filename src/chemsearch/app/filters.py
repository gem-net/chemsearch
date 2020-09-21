from collections import Counter


ALLOWED_FILTERS = tuple()  # set by set_filters_using_config


def get_filters_from_args(args):
    filter_dict = {}
    filter_attrs = [i for i in ALLOWED_FILTERS if i in args]
    for attr in filter_attrs:
        val = args[attr]
        filter_dict[attr] = val
    return filter_dict


def filter_mols(mols, args):
    filter_dict = get_filters_from_args(args)
    if not filter_dict:
        return mols, filter_dict
    for attr, val in filter_dict.items():
        mols = [i for i in mols if getattr(i, attr) == val]
    return mols, filter_dict


def count_filterable(mols):
    counts = dict()
    for attr in ALLOWED_FILTERS:
        attr_counts = Counter([getattr(mol, attr) for mol in mols])
        counts[attr] = attr_counts
    return counts


def update_args(args, attr, val):
    """Get updated copy of args dictionary, popping attr if val is None.

    Aids creation of links that modify URL arguments of current page.
    - add (attr, val) if not present in args.
    - update attr if present in args.
    - if val is None, pop attr.
    """
    new_args = dict(args)
    if 'page' in new_args:
        new_args.pop('page')
    if (attr, val) in args.items():
        new_args.pop(attr)
    else:
        new_args.update({attr: val})
    return new_args


def set_filters_using_config(app):
    global ALLOWED_FILTERS
    if app.config['USE_DRIVE']:
        ALLOWED_FILTERS = (
            'user',
            'category',
        )
    else:
        ALLOWED_FILTERS = (
            'category',
        )
