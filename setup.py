#!/usr/bin/env python
# Copyright (c) 2011-2019, miub developers.
# Distributed under the MIT License. See LICENSE.txt for more info.

"""miubrt - MIUB-radartools - Tools for exploitation of BoXPol radar data and
DWD radar network

"""

import os
import sys
import semver
import warnings
from subprocess import check_output, CalledProcessError

import builtins

DOCLINES = __doc__.split("\n")

CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
Intended Audience :: Education
License :: OSI Approved :: MIT License
Programming Language :: Python
Programming Language :: Python :: 3
Programming Language :: Python :: 3.6
Programming Language :: Python :: 3.7
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3.9
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Atmospheric Science
Operating System :: POSIX :: Linux
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows
"""

NAME = 'miubrt'
MAINTAINER = "Kai Mühlbauer"
MAINTAINER_EMAIL = "kai.muehlbauer@uni-bonn.de"
DESCRIPTION = DOCLINES[0]
LONG_DESCRIPTION = "\n".join(DOCLINES[2:])
URL = "https://git.meteo.uni-bonn.de/projects/miub-radartools"
DOWNLOAD_URL = "https://git.meteo.uni-bonn.de/projects/miub-radartools/repository"
LICENSE = 'MIT'
CLASSIFIERS = filter(None, CLASSIFIERS.split('\n'))
PLATFORMS = ["Linux", "Mac OS-X", "Unix", "Windows"]
MAJOR = 0
MINOR = 1
PATCH = 0
VERSION = '%d.%d.%d' % (MAJOR, MINOR, PATCH)


# Return the git revision as a string
def git_version():
    try:
        git_rev = check_output(['git', 'describe', '--tags', '--long'])
        git_hash = check_output(['git', 'rev-parse', 'HEAD'])
        branch = check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD'])

        git_rev = git_rev.strip().decode('ascii').split('-')
        branch = branch.strip().decode('ascii')

        GIT_REVISION = '-'.join([git_rev[0],
                                 'dev' + git_rev[1]])
        GIT_REVISION = '+'.join([GIT_REVISION,
                                 git_rev[2]])
        GIT_HASH = git_hash.strip().decode('ascii')

        ver = semver.parse_version_info(GIT_REVISION)
        minor = ver.minor
        patch = ver.patch

        if ver.prerelease != 'dev0':
            if 'release' not in branch:
                minor += 1
            else:
                patch += 1
        GIT_REVISION = semver.format_version(ver.major,
                                             minor,
                                             patch,
                                             ver.prerelease)

    except (CalledProcessError, OSError):
        print('MIUBRT: Unable to import git_revision from repository.')
        raise
    return GIT_REVISION, GIT_HASH

# This is a bit hackish: we are setting a global variable so that the main
# miub-radartools __init__ can detect if it is being loaded by the setup
# routine, to # avoid attempting to load components that aren't built yet.
builtins.__MIUBRT_SETUP__ = True


def write_version_py(filename='miubrt/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s
if not release:
    version = full_version
"""
    # Adding the git rev number needs to be done inside write_version_py(),
    # otherwise the import of miubrt.version messes up the build under
    # Python 3.
    SHORT_VERSION = VERSION
    FULL_VERSION = VERSION
    GIT_REVISION = VERSION + '-unknown'
    GIT_HASH = 'unknown'
    ISRELEASED = "'unknown'"

    if os.path.exists('.git'):
        GIT_REVISION, GIT_HASH = git_version()
    elif os.path.exists('miubrt/version.py'):
        # must be a source distribution, use existing version file
        try:
            from miubrt.version import git_revision as GIT_REVISION
            from miubrt.version import git_revision as GIT_HASH
        except ImportError:
            raise ImportError("Unable to import git_revision. Try removing "
                              "miubrt/version.py and the build directory "
                              "before building.")
    else:
        warnings.warn("miubrt source does not contain detailed version info "
                      "via git or version.py, exact version can't be "
                      "retrieved.", UserWarning)
    # parse version using semver
    ver = semver.parse_version_info(GIT_REVISION)

    # get commit count, dev0 means tagged commit -> release
    ISRELEASED = ver.prerelease == 'dev0'
    if not ISRELEASED:
        SHORT_VERSION = semver.format_version(ver.major,
                                              ver.minor,
                                              ver.patch,
                                              ver.prerelease)
        FULL_VERSION = GIT_REVISION

    a = open(filename, 'w')
    try:
        a.write(cnt % {'short_version': SHORT_VERSION,
                       'version': FULL_VERSION,
                       'full_version': GIT_REVISION,
                       'git_revision': GIT_HASH,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()

    return SHORT_VERSION


def setup_package():

    # rewrite version file
    write_version_py()

    from setuptools import setup, find_packages

    with open('requirements.txt', 'r') as f:
        INSTALL_REQUIRES = [rq for rq in f.read().split('\n') if rq != '']

    metadata = dict(
        name=NAME,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        license=LICENSE,
        classifiers=CLASSIFIERS,
        platforms=PLATFORMS,
        install_requires=INSTALL_REQUIRES,
        packages=find_packages(),
        # entry_points={
        #     'console_scripts': [
        #         'radolan2netcdf=miubrt.scripts.radolan2netcdf:main',
        #         'gamic2cfradial=miubrt.scripts.gamic2cfradial:main',
        #         'boxpol2cfradial=miubrt.scripts.boxpol2cfradial:main',
        #         'juxpol2cfradial=miubrt.scripts.juxpol2cfradial:main',
        #         'dwd2cfradial=miubrt.scripts.dwd2cfradial:main',
        #         'create_cvp=miubrt.scripts.create_cvp:main',
        #         'plot_cvp=miubrt.scripts.plot_cvp:main',
        #         'plot_ppi=miubrt.scripts.plot_ppi:main',
        #     ]
        # },
    )

    setup(**metadata)


if __name__ == '__main__':
    setup_package()
