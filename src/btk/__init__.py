"""BlobToolKit."""

import contextlib

from .btk import cli

with contextlib.suppress(ModuleNotFoundError):
    from .lib import pipeline
