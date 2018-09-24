#!/usr/bin/env python3

import pytest
from lib import field

FIELD_ID = 'test_field'
FIELD_DATA = {'values': [1, 2, 3, 4, 5, 4, 3, 2]}
VARIABLE_ID = 'test_variable'
INTEGER_DATA = {'values': [1, 2, 3, 4, 5, 4, 3, 2]}
FLOAT_DATA = {'values': [1.3, 4.2, 6.46733, 4.98, 2.45, 5.4, 3.643, 4.98]}
STRING_DATA = {'values': ['A', 'B', 'third']}
CATEGORY_DATA = {'keys': ['first', 'second'], 'values': [1, 3, 0, 0, 4, 1, 2, 2]}
ARRAY_DATA = {'keys': ['first', 'second'], 'values': [[1, 23, 67, 12], [0, 48, 96, 16]]}


@pytest.fixture(scope='function', name='_my_field')
def my_field():
    return field.Field(FIELD_ID, **FIELD_DATA)


def test_init(_my_field):
    assert _my_field.field_id == FIELD_ID
    assert _my_field.values == FIELD_DATA['values']


def test_get_values_by_index(_my_field):
    my_values = _my_field.get_values_by_indices([2])
    assert isinstance(my_values, list)
    assert my_values == [3]
    assert _my_field.get_values_by_indices(2) == [3]
    assert _my_field.get_values_by_indices([1, 2]) == [2, 3]
    assert _my_field.get_values_by_indices(['string']) == []


def test_get_indices_by_values(_my_field):
    assert _my_field.get_indices_by_values(4) == [3, 5]


def test_get_indices_by_values_with_list(_my_field):
    assert _my_field.get_indices_by_values([3]) == [2, 6]


def test_get_indices_by_values_with_set(_my_field):
    assert _my_field.get_indices_by_values(set([2])) == [1, 7]


# def test_load_yaml(mocker):
#     mocked_read_file = mocker.patch.object(file_io, 'read_file')
#     mocked_read_file.return_value = EXAMPLE_JSON
#     file_io.load_yaml('identifiers.yaml')
#     mocked_read_file.assert_called()
#     mocked_read_file.assert_called_with('identifiers.yaml')
