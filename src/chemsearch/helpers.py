from datetime import timezone

import pandas as pd


def get_utc_naive(dt):
    """Convert timezone aware timestamp to UTC naive timestamp."""
    return dt.astimezone(timezone.utc).replace(tzinfo=None)


def parse_timestamp_str(time):
    """Get naive datetime in UTC."""
    # manual version
    # datetime.strptime(time, '%Y-%m-%dT%H:%M:%S.%fZ').replace(tzinfo=timezone.utc)
    timestamp = pd.to_datetime(time).replace(tzinfo=timezone.utc)
    # dt = timestamp.to_pydatetime().replace(tzinfo=timezone.utc)
    return get_utc_naive(timestamp)


def to_google_time(timestamp):
    return timestamp.strftime('%Y-%m-%dT%H:%M:%S.%fZ')
