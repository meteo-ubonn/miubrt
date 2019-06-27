#!/usr/bin/env python
# Copyright (c) 2019-2020, miub developers.
# Distributed under the MIT License. See LICENSE.txt for more info.

"""
miub-radartools
===============

"""

# Detect if we're being called as part of setup procedure
try:
    __MIUBRT_SETUP__
except NameError:
    __MIUBRT_SETUP__ = False

if __MIUBRT_SETUP__:
    import sys as _sys
    _sys.stderr.write("Running from source directory.\n")
    del _sys
else:
    # Make sure that deprecation warnings get printed by default
    import warnings as _warnings
    _warnings.filterwarnings("always", category=DeprecationWarning,
                             module='miubrt')

    # versioning
    from .version import git_revision as __git_revision__  # noqa
    from .version import version as __version__  # noqa

    # packages
    from . import util  # noqa
    from . import io  # noqa
    from . import xarray  # noqa