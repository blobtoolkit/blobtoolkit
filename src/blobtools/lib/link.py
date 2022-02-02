#!/usr/bin/env python3
"""Add links to metadata."""

import re
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import urlopen


def _expand_link(url, values, identifiers):
    """Expand a link url using dict of values."""
    parts = re.split(r'\{|\}', url)
    for i in range(1, len(parts), 2):
        if parts[i] == 'id':
            parts[i] = identifiers[0]
        else:
            parts[i] = str(values[parts[i]])
    return ''.join(parts)


def _test_link(url):
    """Check url exists."""
    url = quote(url, ':/')
    try:
        urlopen(url)
    except HTTPError as error:
        print('HTTPError: {}'.format(error.code))
    except URLError as error:
        print('URLError: {}'.format(error.reason))
    else:
        return True


def add(string, meta, identifiers, skip_test):
    """Add a link to meta."""
    path, url = string.split('=', 1)
    keys = path.split('.')
    links = meta.links
    values = meta
    if keys[0] == 'position':
        skip_test = True
    if skip_test:
        valid = True
    else:
        for key in keys[:-2]:
            if isinstance(values, dict):
                values = values[key]
            else:
                values = getattr(values, key)
        valid = _test_link(_expand_link(url, values, identifiers))
    if valid:
        for key in keys[:-1]:
            if key not in links:
                links[key] = {}
            links = links[key]
        links.update({keys[-1]: url})
