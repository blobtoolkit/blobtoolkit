#!/usr/bin/env python3
"""
Field Class module.

Generic, identifier, variable, category, array and object field types.
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
        self.update_data(**kwargs)
        # for key, value in kwargs.items():
        #     setattr(self, key, value)

    def update_data(self, **kwargs):
        """Update values and keys for an existing field."""
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
                values = [values]
            else:
                values = set(values)
        print(values)
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


class Identifier(Field):
    """Class for record identifiers."""

    def __init__(self, field_id, **kwargs):
        """Init Identifier class."""
        super().__init__(field_id, **kwargs)
        self.type = 'identifier'

    @staticmethod
    def check_unique(entries):
        """Check all entries are unique."""
        unique = set()
        for entry in entries:
            unique.add(entry)
        return len(unique) == len(entries)

    def validate_list(self, names):
        """Ensure all list entries are unique and match identifiers."""
        valid = True
        if self.check_unique(names):
            if len(self.get_indices_by_values(names)) != len(names):
                valid = False
        else:
            valid = False
        return valid


class Variable(Field):
    """Class for variable field type."""

    def __init__(self, field_id, **kwargs):
        """Init Variable class."""
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


__all__ = ['Field', 'Identifier', 'Variable']


if __name__ == '__main__':
    import doctest
    doctest.testmod()
