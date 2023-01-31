#!/usr/bin/env python
# Copyright (c) 2019-2023, miub developers.
# Distributed under the MIT License. See LICENSE.txt for more info.

"""
miub-radartools
===============

"""

# Make sure that deprecation warnings get printed by default
import warnings as _warnings

_warnings.filterwarnings("always", category=DeprecationWarning, module="miubrt")

# versioning
try:
    from .version import version as __version__
except Exception:
    # Local copy or not installed with setuptools.
    # Disable minimum version checks on downstream libraries.
    __version__ = "999"

# packages
from . import io  # noqa
from . import ipol  # noqa
from . import plot  # noqa
from . import util  # noqa
from . import xarray  # noqa
