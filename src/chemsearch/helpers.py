from datetime import timezone

import bleach
import pandas as pd


def get_utc_naive(dt):
    """Convert timezone aware timestamp to UTC naive timestamp."""
    return dt.astimezone(timezone.utc).replace(tzinfo=None)


def parse_timestamp_str(time):
    """Get naive datetime in UTC."""
    # manual version
    # datetime.strptime(time, '%Y-%mol-%dT%H:%M:%S.%fZ').replace(tzinfo=timezone.utc)
    timestamp = pd.to_datetime(time).replace(tzinfo=timezone.utc)
    # dt = timestamp.to_pydatetime().replace(tzinfo=timezone.utc)
    return get_utc_naive(timestamp)


def to_google_time(timestamp):
    return timestamp.strftime('%Y-%mol-%dT%H:%M:%S.%fZ')


def clean_html(html):
    allowed_tags = ['a', 'abbr', 'acronym', 'b', 'blockquote', 'code',
                    'em', 'i', 'li', 'ol', 'strong', 'ul', 'img',
                    'pre', 'h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'p', 'hr',
                    'table', 'thead', 'tr', 'th', 'tbody', 'td']
    allowed_attrs = {'ol': ['start']}
    html_clean = bleach.linkify(bleach.clean(html, tags=allowed_tags,
                                             attributes=allowed_attrs))
    html_clean = html_clean.replace(
        '<table>', '<table class="table table-responsive table-hover">')
    return html_clean
