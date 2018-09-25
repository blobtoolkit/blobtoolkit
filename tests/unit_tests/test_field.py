#!/usr/bin/env python3

import pytest
from lib import Field, Variable

FIELD_ID = 'test_field'
FIELD_DATA = {'keys': ['first', 'second'], 'values': [1, 2, 3, 4, 5, 4, 3, 2]}
VARIABLE_ID = 'test_variable'
INTEGER_DATA = {'values': [1, 2, 3, 4, 5, 4, 3, 2]}
FLOAT_DATA = {'values': [1.3, 4.2, 6.46733, 4.98, 2.45, 5.4, 3.643, 4.98]}
STRING_DATA = {'values': ['A', 'B', 'third']}
CATEGORY_DATA = {'keys': ['first', 'second', 'third', 'fourth', 'fifth'],
                 'values': [1, 3, 0, 0, 4, 1, 2, 2]}
ARRAY_DATA = {'keys': ['first', 'second'], 'values': [[1, 23, 67, 12], [0, 48, 96, 16]]}


@pytest.fixture(scope='function', name='_my_field')
def my_field():
    return Field(FIELD_ID, **FIELD_DATA)


def test_field_init(_my_field):
    assert _my_field.field_id == FIELD_ID
    assert _my_field.values == FIELD_DATA['values']
    assert _my_field.keys == FIELD_DATA['keys']
    assert _my_field.type == 'generic'


def test_get_values_by_index(_my_field):
    my_values = _my_field.get_values_by_indices([2])
    assert isinstance(my_values, list)
    assert my_values == [3]
    assert _my_field.get_values_by_indices(2) == [3]
    assert _my_field.get_values_by_indices([1, 2]) == [2, 3]
    assert _my_field.get_values_by_indices(['string']) == []


def test_get_indices_by_values(_my_field):
    my_indices = _my_field.get_indices_by_values(4)
    assert isinstance(my_indices, list)
    assert my_indices == [3, 5]
    assert _my_field.get_indices_by_values([3]) == [2, 6]
    assert _my_field.get_indices_by_values(set([2])) == [1, 7]


def test_subset(_my_field):
    assert _my_field.subset == FIELD_DATA['values']
    _my_field.subset = [1, 3, 5]
    assert _my_field.subset == [1, 3, 5]
    _my_field.subset = []
    assert _my_field.subset == []
    _my_field.subset = False
    assert _my_field.subset == FIELD_DATA['values']


def test_select_records(_my_field, mocker):
    mocked_get_values_by_indices = mocker.patch.object(Field, 'get_values_by_indices')
    mocked_get_values_by_indices.return_value = [1, 3, 5]
    values = _my_field.select_records([0, 2, 4])
    mocked_get_values_by_indices.assert_called()
    mocked_get_values_by_indices.assert_called_with([0, 2, 4])
    assert values == [1, 3, 5]
    assert _my_field.subset == [1, 3, 5]
    mocked_get_values_by_indices.return_value = []
    values = _my_field.select_records([])
    assert values == []
    values = _my_field.select_records(False)
    assert values is False


@pytest.fixture(scope='function', name='_my_variable')
def my_variable():
    return Variable(VARIABLE_ID, **INTEGER_DATA)


def test_variable_init(_my_variable):
    assert _my_variable.field_id == VARIABLE_ID
    assert _my_variable.values == INTEGER_DATA['values']
    assert _my_variable.type == 'variable'


def test_get_indices_in_range(_my_variable):
    assert _my_variable.get_indices_in_range(2, 3) == -1
    assert _my_variable.get_indices_in_range([2]) == -2
    assert _my_variable.get_indices_in_range(['1', '4']) == -3
    my_indices = _my_variable.get_indices_in_range([2, 4])
    assert isinstance(my_indices, list)
    assert my_indices == [1, 2, 3, 5, 6, 7]
    assert _my_variable.get_indices_in_range([2.3, 3.4]) == [2, 6]
    assert _my_variable.get_indices_in_range([2.3, 2.5]) == []
    assert _my_variable.get_indices_in_range([2, 4], True) == [0, 4]


def test_sum_values(_my_variable):
    assert _my_variable.sum_values() == 24
    _my_variable.subset = [1, 3, 5]
    assert _my_variable.sum_values() == 9
    _my_variable.subset = False
    assert _my_variable.sum_values() == 24
    _my_variable.subset = []
    assert _my_variable.sum_values() == 0


# def test_load_yaml(mocker):
#     mocked_read_file = mocker.patch.object(file_io, 'read_file')
#     mocked_read_file.return_value = EXAMPLE_JSON
#     file_io.load_yaml('identifiers.yaml')
#     mocked_read_file.assert_called()
#     mocked_read_file.assert_called_with('identifiers.yaml')
