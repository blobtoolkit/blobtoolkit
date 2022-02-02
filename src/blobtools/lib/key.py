#!/usr/bin/env python3
"""Add key values to metadata."""

import re


def add(string, meta, replace):
    """Add a key value to meta."""
    path, value = string.split('=')
    if re.match(r'^\[.+\]$', value):
        value = value.strip('[]').split(',')
    keys = path.split('.')
    if len(keys) == 1:
        setattr(meta, keys[0], value)
        return True
    if hasattr(meta, keys[0]):
        current = getattr(meta, keys[0])
    else:
        current = {}
        setattr(meta, keys[0], current)
    for key in keys[1:-1]:
        if key not in current:
            current[key] = {}
        current = current[key]
    if not replace and keys[-1] in current:
        if isinstance(current[keys[-1]], list):
            if isinstance(value, list):
                value = current[keys[-1]] + value
            else:
                val = value
                value = current[keys[-1]][:]
                value.append(val)
        elif current[keys[-1]]:
            value = [current[keys[-1]], value]
    current.update({keys[-1]: value})
    return True
