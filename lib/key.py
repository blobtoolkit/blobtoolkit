#!/usr/bin/env python3
"""Add key values to metadata."""


def add(string, meta):
    """Add a key value to meta."""
    path, value = string.split('=')
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
    current.update({keys[-1]: value})
    return True
