#!/usr/bin/env python3
# pylint: disable=wrong-import-position
import pytest

from lib import Metadata, Taxdump, taxid

TAXON_ID = 3702
DATASET_ID = 'DS1'
DATASET_DATA = {'taxon': {'common_name': 'thale-cress'}}
TAXDUMP_DIR = '/dev/null'
TAXDUMP_DATA = {'names': {}, 'ancestors': {}, 'ranks': {}}


@pytest.fixture(scope='function', name='_my_taxdump')
def my_taxdump():
    return Taxdump(TAXDUMP_DIR, **TAXDUMP_DATA)


@pytest.fixture(scope='function', name='_my_meta')
def my_meta():
    return Metadata(DATASET_ID, **DATASET_DATA)


# def mock_lineage(mocker):
#     return mocked_lineage
#
#
# def mock_no_lineage(mocker):
#     mocked_lineage = mocker.patch.object(Taxdump, 'lineage')
#     mocked_lineage.return_value = {}
#     return mocked_lineage


def test_add_invalid(_my_taxdump, _my_meta, mocker):
    mocked_lineage = mocker.patch.object(Taxdump, 'lineage')
    mocked_lineage.return_value = {}
    taxon = taxid.add(0, _my_taxdump, _my_meta)
    assert isinstance(taxon, dict)
    assert 'taxid' not in taxon
    assert taxon['common_name'] == DATASET_DATA['taxon']['common_name']


def test_add(_my_taxdump, _my_meta, mocker):
    mocked_lineage = mocker.patch.object(Taxdump, 'lineage')
    mocked_lineage.return_value = {'superkingdom': 'Eukaryota',
                                   'genus': 'Arabidopsis',
                                   'species': 'Arabidopsis thaliana'}
    taxon = taxid.add(TAXON_ID, _my_taxdump, _my_meta)
    assert isinstance(taxon, dict)
    assert taxon['taxid'] == TAXON_ID
    assert taxon['superkingdom'] == DATASET_DATA['taxon']['superkingdom']
    assert taxon['common_name'] == DATASET_DATA['taxon']['common_name']
