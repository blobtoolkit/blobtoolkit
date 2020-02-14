"""Configuration parameters for pytest."""
# pylint: disable=unused-import
import sys
from os.path import abspath
from os.path import dirname as d
from os.path import join

from pytest_mock import mocker

ROOT_DIR = d(d(abspath(__file__)))
sys.path.insert(0, ROOT_DIR)
sys.path.insert(0, abspath(join(d(__file__), 'lib')))
print(sys.path)
