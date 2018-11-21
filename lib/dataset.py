#!/usr/bin/env python3
"""Dataset Class module."""


class Metadata():
    """Class for field and dataset metadata."""

    __slots__ = ['dataset_id',
                 'assembly',
                 '_field_ids',
                 '_field_list',
                 'fields',
                 'links',
                 'name',
                 'plot',
                 'reads',
                 'record_type',
                 'records',
                 'settings',
                 'similarity',
                 'static_fields',
                 'taxon'
                 ]

    def __init__(self, dataset_id, **kwargs):
        """Init Dataset class."""
        self.dataset_id = dataset_id
        self.fields = []
        self.name = ''
        self.links = {}
        self.record_type = 'record'
        self.update_data(**kwargs)
        self._field_list = self.list_fields()
        self._field_ids = list(self._field_list.keys())

    def update_data(self, **kwargs):
        """Update values and keys for an existing field."""
        for key, value in kwargs.items():
            if key == 'id':
                key = 'dataset_id'
            setattr(self, key, value)

    def list_fields(self):
        """List all fields in dataset."""
        field_ids = self._list_fields(self.fields)
        return field_ids

    def has_field(self, field_id):
        """Return true if a field_id is present."""
        return field_id in self._field_ids

    def add_parents(self, parents):
        """Add field metadata."""
        parent = self._field_list
        fields = self.fields
        for required in parents:
            if isinstance(required, dict):
                if isinstance(parent, dict):
                    if required['id'] not in parent:
                        parent.update({required['id']: required})
                        fields.append(required)
                    parent = parent[required['id']]
                elif required['id'] in self._field_list:
                    parent = self._field_list[required['id']]
                else:
                    parent.append(required)
                    parent = required
            else:
                if required not in parent:
                    parent.update({required: []})
                parent = parent[required]
        return parent

    def add_field(self, parents=None, **kwargs):
        """Add field metadata."""
        if parents is None:
            parents = []
        parent = self.add_parents(parents)
        if self.has_field(kwargs['field_id']):
            meta = self._field_list[kwargs['field_id']]
        else:
            meta = {}
            parent.append(meta)
        for key, value in kwargs.items():
            if key == 'field_id':
                key = 'id'
            meta[key] = value

    def to_dict(self):
        """Create a dict of metadata."""
        data = {}
        for key in self.__slots__:
            if not key.startswith('_'):
                if hasattr(self, key):
                    data[key] = getattr(self, key)
        return data

    @staticmethod
    def _list_fields(parent, fields=None):
        if fields is None:
            fields = {}
        for field in parent:
            fields.update({field['id']: field})
            if 'children' in field:
                Metadata._list_fields(field['children'], fields)
            if 'data' in field:
                Metadata._list_fields(field['data'], fields)
        return fields


if __name__ == '__main__':
    import doctest
    doctest.testmod()
