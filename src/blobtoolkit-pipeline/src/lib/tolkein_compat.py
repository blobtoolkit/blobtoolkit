"""Compatibility imports for tolkein modules."""

import warnings

try:
    from genomehubs.tolkein import tofetch
    from genomehubs.tolkein import tofile
    from genomehubs.tolkein import tolog
except ModuleNotFoundError:
    try:
        from genomehubs.vendor.tolkein import tofetch
        from genomehubs.vendor.tolkein import tofile
        from genomehubs.vendor.tolkein import tolog
    except ModuleNotFoundError:
        from tolkein import tofetch
        from tolkein import tofile
        from tolkein import tolog

        warnings.warn(
            "Importing from `tolkein` is deprecated. Please install `genomehubs` "
            "and use `genomehubs.tolkein`.",
            DeprecationWarning,
            stacklevel=2,
        )

__all__ = ["tofetch", "tofile", "tolog"]
