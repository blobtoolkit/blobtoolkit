#!/usr/bin/env python3
"""
Example Class module.

Tests are included in docstrings and run using doctest.
"""


class Field():
    """Parent class for specific field types."""

    __slots__ = ['field_id', 'values', 'keys']

    def __init__(self, field_id, **kwargs):
        """Init Field class."""
        self.field_id = field_id
        self.values = []
        self.keys = []
        for key, value in kwargs.items():
            setattr(self, key, value)

    def get_values_by_indices(self, indices):
        """Get values for all records at matching indices."""
        values = []
        if not isinstance(indices, list):
            indices = [indices]
        for index in indices:
            try:
                values.append(self.values[index])
            except TypeError:
                pass
        return values

    def get_indices_by_value(self, values):
        """Get indices for all records matching a value."""
        if not isinstance(values, set):
            if not isinstance(values, list):
                values = set([values])
            else:
                values = set(values)
        indices = [i for i, x in enumerate(self.values) if x in values]
        return indices


if __name__ == '__main__':
    import doctest
    doctest.testmod()
