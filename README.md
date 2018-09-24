## Add new datatypes to Blobtools

[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![Build Status](https://travis-ci.org/blobtoolkit/blobtools-add.svg?branch=master)](https://travis-ci.org/blobtoolkit/blobtools-add)
[![Coverage Status](https://coveralls.io/repos/github/blobtoolkit/blobtools-add/badge.svg?branch=master)](https://coveralls.io/github/blobtoolkit/blobtools-add?branch=master)

Repository based on [CI template](http://github.com/rjchallis/test)

### About

This repository is a space to develop new functionality to bridge the gap between [Blobtools](https://github.com/DRL/blobtools) and the [BlobToolKit Viewer](https://github.com/blobtoolkit/viewer).

### Installing

Clone/Fork the repository:
```
git clone https://github.com/blobtoolkit/blobtools-add
```

Install dependencies:
```
cd test
pip install -r requirements.txt
```

### Contributing

If you find a problem or want to suggest a feature, please submit an issue.

If you want to contribute code, pull requests are welcome but please make sure your code passes the linting, testing and style checks before submitting a pull request.

Run linter/testing locally:
```
./run_tests.sh
```

Set up pre-commit hook to automate test running:
```
ln -s ../../pre-commit.sh .git/hooks/pre-commit
```
