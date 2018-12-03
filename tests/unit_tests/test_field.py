#!/usr/bin/env python3
# pylint: disable=wrong-import-position
import pytest
from lib import Array, Category, Field, Identifier, Variable

FIELD_ID = 'test_field'
FIELD_DATA = {'keys': ['first', 'second'], 'values': [1, 2, 3, 4, 5, 4, 3, 2]}
STRING_DATA = {'values': ['A', 'B', 'third']}


@pytest.fixture(scope='function', name='_my_field')
def my_field():
    return Field(FIELD_ID, **FIELD_DATA)


def test_update_values(_my_field):
    _my_field.update_values([10, 20, 30, 40, 50, 40, 30, 20])
    assert _my_field.values[2] == 30
    _my_field.update_values([110, 120, 130, 140, 150, 140])
    assert _my_field.values[2] == 30


def test_update_keys(_my_field):
    _my_field.update_keys(['key1', 'key2'])
    assert _my_field.keys[1] == 'key2'


def test_field_init(_my_field):
    assert _my_field.field_id == FIELD_ID
    assert _my_field.values == FIELD_DATA['values']
    assert _my_field.keys == FIELD_DATA['keys']
    assert _my_field.type == 'field'


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


def test_update_data(_my_field):
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
CATEGORY_ID = 'test_category'
CATEGORY_DATA = {'keys': ['first', 'second', 'third', 'fourth', 'fifth'],
                 'values': [1, 3, 0, 0, 4, 1, 2, 2]}
CATEGORY_VALUES_ONLY = {'values': ['A', 'C', 'B', 'D', 'A', 'D', 'B', 'E']}
# TODO: convert list to category


@pytest.fixture(scope='function', name='_my_category')
def my_category():
    return Category(CATEGORY_ID, **CATEGORY_DATA)


def test_category_init(_my_category):
    assert _my_category.field_id == CATEGORY_ID
    assert _my_category.values == CATEGORY_DATA['values']
    assert _my_category.type == 'category'


def test_category_init_without_keys():
    keyless_category = Category(CATEGORY_ID, **CATEGORY_VALUES_ONLY)
    assert keyless_category.values == [0, 1, 2, 3, 0, 3, 2, 4]
    assert keyless_category.keys == ['A', 'C', 'B', 'D', 'E']


ARRAY_ID = 'test_array'
ARRAY_DATA = (
    {
        'keys': ['first', 'second'],
        'values': [
            [1, 23, 67, 12],
            [0, 48, 96, 16],
            [1, 45, 72, 18],
            [0, 41, 88, 12]
        ]
    }
)
ARRAY_VALUES_ONLY = (
    {
        'category_slot': 0,
        'values': [
            ['second', 23, 67, 12],
            ['first', 48, 96, 16],
            ['second', 45, 72, 18],
            ['first', 41, 88, 12]
        ]
    }
)


@pytest.fixture(scope='function', name='_my_array')
def my_array():
    return Array(ARRAY_ID, **ARRAY_DATA)


def test_array_init(_my_array):
    assert _my_array.field_id == ARRAY_ID
    assert _my_array.values == ARRAY_DATA['values']
    assert _my_array.type == 'array'


def test_array_init_without_keys():
    keyless_array = Array(ARRAY_ID, **ARRAY_VALUES_ONLY)
    assert keyless_array.values[3][0] == 1
    assert keyless_array.keys == ['second', 'first']


def test_get_values_by_indices_for_slots(_my_array):
    values = _my_array.get_values_by_indices_for_slots([2, 3], [1, 2])
    assert isinstance(values, list)
    assert len(values) == 2
    assert isinstance(values[0], list)
    assert values[0] == [45, 72]
    values = _my_array.get_values_by_indices_for_slots([0, 1], 3)
    assert len(values) == 2
    assert values == [12, 16]


def test_update_slots(_my_array):
    _my_array.update_slots([0, 1, 0, 1])
    assert _my_array.values[2][0] == 0
    _my_array.update_slots([2, 4, 6, 7], slot=3)
    assert _my_array.values[1][3] == 4


# TODO: test Object
