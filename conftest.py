import sys
from pytest_mock import mocker
from os.path import dirname as d
from os.path import abspath, join
root_dir = d(d(abspath(__file__)))
sys.path.insert(0, root_dir)
sys.path.insert(0, abspath(join(d(__file__), 'lib')))
