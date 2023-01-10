#!/usr/bin/env python3
"""
Fetch BlobToolKit Pipeline data.

NOT YET IMPLEMENTED

Usage: blobtoolkit-pipeline data --config YAML

Options:
    --config YAML  YAML format configuration filename.
"""

import logging
import re
import sys
from collections import defaultdict

from docopt import DocoptExit
from docopt import docopt

logger_config = {
    "level": logging.INFO,
    "format": "%(asctime)s [%(levelname)s] line %(lineno)d %(message)s",
    "filemode": "w",
}
logging.basicConfig(**logger_config)
logger = logging.getLogger()


def main(rename=None):
    """Entry point."""
    docs = __doc__
    if rename is not None:
        docs = docs.replace("blobtoolkit-pipeline", rename)
    try:
        args = docopt(docs)
    except DocoptExit as e:
        raise DocoptExit from e
