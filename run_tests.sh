#!/bin/sh

# lint code in lib directory
pylint --rcfile=.pylintrc lib -f parseable -r n &&
# check codestyle
pycodestyle lib --max-line-length=120 &&
# check docstyle
pydocstyle lib &&

# lint code in test directory with alternate config
pylint --rcfile=.pylintrc_tests tests/unit_tests/* -f parseable -r n &&
# check codestyle for tests
pycodestyle tests/unit_tests/* --max-line-length=120 &&
# requiring docstrings for tests would be excessive so don't run pydocstyle

# run tests and generate coverage report
py.test --cov-config .coveragerc --doctest-modules --cov=lib --cov-report term-missing
