#!/usr/bin/env python3

import pytest
from lib import field

FIELD_ID = 'test_field'
FIELD_DATA = {'values': [1, 2, 3, 4, 5, 4, 3, 2]}


@pytest.fixture(scope='function', name='_my_field')
def my_field():
    return field.Field(FIELD_ID, **FIELD_DATA)


def test_init(_my_field):
    assert _my_field.field_id == FIELD_ID
    assert _my_field.values == FIELD_DATA['values']


def test_get_value_by_index(_my_field):
    assert _my_field.get_value_by_index(2) == 3


def test_get_indices_by_value(_my_field):
    assert _my_field.get_indices_by_value(4) == [3, 5]


def test_get_indices_by_value_with_list(_my_field):
    assert _my_field.get_indices_by_value([3]) == [2, 6]


def test_get_indices_by_value_with_set(_my_field):
    assert _my_field.get_indices_by_value(set([2])) == [1, 7]


# def test_load_yaml(mocker):
#     mocked_read_file = mocker.patch.object(file_io, 'read_file')
#     mocked_read_file.return_value = EXAMPLE_JSON
#     file_io.load_yaml('identifiers.yaml')
#     mocked_read_file.assert_called()
#     mocked_read_file.assert_called_with('identifiers.yaml')
