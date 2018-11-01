#!/usr/bin/env python3

import pytest
from lib import Field, Identifier, Variable

FIELD_ID = 'test_field'
FIELD_DATA = {'keys': ['first', 'second'], 'values': [1, 2, 3, 4, 5, 4, 3, 2]}
STRING_DATA = {'values': ['A', 'B', 'third']}


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


def test_update_values(_my_field):
    _my_field.update_data(**STRING_DATA)
    assert _my_field.values == STRING_DATA['values']
    assert _my_field.get_indices_by_values('B') == [1]
    assert _my_field.get_indices_by_values([4]) == []


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


IDENTIFIER_ID = 'identifiers'
IDENTIFIERS = {'values': ['r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8']}
DUPLICATED_IDENTIFIERS = {'values': ['r1', 'r2', 'r3', 'r3', 'r5', 'r6', 'r7', 'r8']}


@pytest.fixture(scope='function', name='_my_identifier')
def my_identifier():
    return Identifier(IDENTIFIER_ID, **IDENTIFIERS)


def test_identifier_init(_my_identifier):
    assert _my_identifier.field_id == IDENTIFIER_ID
    assert _my_identifier.values == IDENTIFIERS['values']
    assert _my_identifier.type == 'identifier'


def test_identifier_check_unique(_my_identifier):
    assert _my_identifier.check_unique(['a', 'b', 'c']) is True
    assert _my_identifier.check_unique(['a', 'b', 'b']) is False


def test_identifier_validate_list(_my_identifier):
    assert _my_identifier.validate_list(['r1', 'r2', 'r3']) is True
    assert _my_identifier.validate_list(['r1', 'r2', 'r2']) is False
    assert _my_identifier.validate_list(['r1', 'r2', 's3']) is False


VARIABLE_ID = 'test_variable'
INTEGER_DATA = {'values': [1, 2, 3, 4, 5, 4, 3, 2]}
# TODO: test Variable with float data
# FLOAT_DATA = {'values': [1.3, 4.2, 6.46733, 4.98, 2.45, 5.4, 3.643, 4.98]}


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


def test_sum_values(_my_variable, mocker):
    assert _my_variable.sum_values() == 24
    mocked_prop = mocker.PropertyMock(return_value=[1, 3, 5])
    type(_my_variable).subset = mocked_prop
    assert _my_variable.sum_values() == 9
    mocked_prop = mocker.PropertyMock(return_value=[])
    type(_my_variable).subset = mocked_prop
    assert _my_variable.sum_values() == 0


# TODO: test Category
CATEGORY_DATA = {'keys': ['first', 'second', 'third', 'fourth', 'fifth'],
                 'values': [1, 3, 0, 0, 4, 1, 2, 2]}
# TODO: convert list to category

# TODO: test Array
ARRAY_DATA = {'keys': ['first', 'second'], 'values': [[1, 23, 67, 12], [0, 48, 96, 16]]}

# TODO: test Object
