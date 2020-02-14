#!/usr/bin/env python3
# pylint: disable=wrong-import-position
import pytest
from lib import Taxdump

TAXON_ID = 3888
DATASET_ID = 'DS2'
DATASET_DATA = {'taxon': {'common_name': 'pea'}}
TAXDUMP_DIR = '/dev/null'
TAXDUMP_DATA = {'ancestors': {}, 'names': {}, 'ranks': {}}


@pytest.fixture(scope='function', name='_my_taxdump')
def my_taxdump():
    return Taxdump(TAXDUMP_DIR, **TAXDUMP_DATA)


# def test_add(_my_taxdump, _my_meta, mocker):
#     mocked_lineage = mocker.patch.object(Taxdump, 'lineage')
#     mocked_lineage.return_value = {'superkingdom': 'Eukaryota',
#                                    'genus': 'Arabidopsis',
#                                    'species': 'Arabidopsis thaliana'}
#     taxon = taxid.add(TAXON_ID, _my_taxdump, _my_meta)
#     assert isinstance(taxon, dict)
#     assert taxon['taxid'] == TAXON_ID
#     assert taxon['superkingdom'] == DATASET_DATA['taxon']['superkingdom']
#     assert taxon['common_name'] == DATASET_DATA['taxon']['common_name']
