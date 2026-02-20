"""BlobToolKit."""

import contextlib

from .blobtools import cli
from .lib import add
from .lib import filter
from .lib.version import __version__

with contextlib.suppress(ModuleNotFoundError):
    from .lib import host
    from .lib import view

# from .field import Array, Category, Field, Identifier, MultiArray, Variable
# from .taxdump import Taxdump
# from .dataset import Metadata
