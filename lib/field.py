#!/usr/bin/env python3
"""
Example Class module.

Tests are included in docstrings and run using doctest.
"""

from functools import reduce
import numbers
from operator import add


class Field():
    """Parent class for specific field types."""

    __slots__ = ['field_id', 'values', 'keys', 'type', '_subset']

    def __init__(self, field_id, **kwargs):
        """Init Field class."""
        self.field_id = field_id
        self.type = 'generic'
        self.values = []
        self.keys = []
        self._subset = False
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

    def get_indices_by_values(self, values):
        """Get indices for all records matching a value."""
        if not isinstance(values, set):
            if not isinstance(values, list):
                values = set([values])
            else:
                values = set(values)
        indices = [i for i, x in enumerate(self.values) if x in values]
        return indices

    def select_records(self, indices=False):
        """Create a subset of values from a list of indices."""
        values = []
        if not isinstance(indices, list):
            values = False
        else:
            values = self.get_values_by_indices(indices)
        self.subset = values
        return values

    @property
    def subset(self):
        """Subset property."""
        subset = []
        if isinstance(self._subset, list):
            subset = self._subset
        else:
            subset = self.values
        return subset

    @subset.setter
    def subset(self, values=False):
        """Set subset using indices."""
        if isinstance(values, list):
            self._subset = values
        else:
            self._subset = False


class Variable(Field):
    """Class for variable field type."""

    def __init__(self, field_id, **kwargs):
        """Init Field class."""
        super().__init__(field_id, **kwargs)
        self.type = 'variable'

    def get_indices_in_range(self, min_max, invert=False):
        """Get indices for all records in range."""
        if not isinstance(min_max, list):
            return -1
        if not len(min_max) == 2:
            return -2
        if not all(isinstance(x, numbers.Real) for x in min_max):
            return -3
        indices = []
        for index, value in enumerate(self.values):
            if invert:
                if value < min_max[0] or value > min_max[1]:
                    indices.append(index)
            else:
                if min_max[0] <= value <= min_max[1]:
                    indices.append(index)
        return indices

    def sum_values(self):
        """Sum values in subset."""
        try:
            total = reduce(add, self.subset)
        except TypeError:
            total = 0
        return total


__all__ = ['Field', 'Variable']


if __name__ == '__main__':
    import doctest
    doctest.testmod()
