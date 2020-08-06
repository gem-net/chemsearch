import math

from flask import current_app, abort
from more_itertools import chunked


def get_page_count(n_items: int):
    per_page = current_app.config['MOLECULES_PER_PAGE']
    return math.ceil(n_items / per_page)


def get_page_items_or_404(iterable, page_no):
    if page_no == 1 and not iterable:
        return iterable
    items = []
    try:
        per_page = current_app.config['MOLECULES_PER_PAGE']
        items = list(chunked(iterable, per_page))[page_no - 1]
    except IndexError:
        abort(404)
    return items
