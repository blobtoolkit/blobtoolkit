#!/usr/bin/env python3

from os.path import abspath
import binascii
from lib import file_io


VALUE = ['record1', 'record2', 'record3', 'record4', 'record5']
EXAMPLE_JSON = '[\n "record1",\n "record2",\n "record3",\n "record4",\n "record5"\n]'
EXAMPLE_YAML = '- record1\n- record2\n- record3\n- record4\n- record5\n'


def test_read_file_that_exists():
    assert file_io.read_file('tests/files/infile') == 'testfile content\n'
    assert file_io.read_file('tests/files/infile.gz') == 'testfile content\n'


def test_read_file_that_does_not_exist():
    assert file_io.read_file('nofile') is None
    assert file_io.read_file('nofile.gz') is None


def test_load_yaml(mocker):
    mocked_read_file = mocker.patch.object(file_io, 'read_file')
    mocked_read_file.return_value = EXAMPLE_JSON
    file_io.load_yaml('identifiers.yaml')
    mocked_read_file.assert_called()
    mocked_read_file.assert_called_with('identifiers.yaml')


def test_write_file_json(tmpdir):
    path = tmpdir.mkdir('sub').join('outfile.json')
    pathname = abspath(path)
    file_io.write_file(pathname, VALUE)
    assert len(tmpdir.listdir()) == 1
    assert path.read() == EXAMPLE_JSON
    tmpdir.remove()


def test_write_file_yaml(tmpdir):
    path = tmpdir.mkdir('sub').join('outfile.yaml')
    pathname = abspath(path)
    file_io.write_file(pathname, VALUE)
    assert len(tmpdir.listdir()) == 1
    assert path.read() == EXAMPLE_YAML
    tmpdir.remove()


def test_write_file_plain(tmpdir):
    path = tmpdir.mkdir('sub').join('outfile.txt')
    pathname = abspath(path)
    example_string = 'file content\n'
    file_io.write_file(pathname, example_string)
    assert len(tmpdir.listdir()) == 1
    assert path.read() == example_string
    tmpdir.remove()


def test_write_file_invalid_path():
    example_string = 'file content\n'
    assert file_io.write_file('path/does/not/exist', example_string) is False
    assert file_io.write_file('path/does/not/exist.gz', example_string) is False


def test_write_file_json_gz(tmpdir):
    path = tmpdir.mkdir('sub').join('outfile.json.gz')
    pathname = abspath(path)
    file_io.write_file(pathname, VALUE)
    assert len(tmpdir.listdir()) == 1
    with open(pathname, 'rb') as fh:
        assert binascii.hexlify(fh.read(2)) == b'1f8b'
    tmpdir.remove()


def test_write_file_yaml_gz(tmpdir):
    path = tmpdir.mkdir('sub').join('outfile.yaml.gz')
    pathname = abspath(path)
    file_io.write_file(pathname, VALUE)
    assert len(tmpdir.listdir()) == 1
    with open(pathname, 'rb') as fh:
        assert binascii.hexlify(fh.read(2)) == b'1f8b'
    tmpdir.remove()
