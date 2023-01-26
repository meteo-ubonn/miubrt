#!/usr/bin/env python
# Copyright (c) 2023, miub developers.
# Distributed under the MIT License. See LICENSE.txt for more info.

"""miubrt - MIUB-radartools - Tools for exploitation of BoXPol radar data and
DWD radar network

"""

DOCLINES = __doc__.split("\n")

CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
Intended Audience :: Education
License :: OSI Approved :: MIT License
Programming Language :: Python
Programming Language :: Python :: 3
Programming Language :: Python :: 3.9
Programming Language :: Python :: 3.10
Programming Language :: Python :: 3.11
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Atmospheric Science
Operating System :: POSIX :: Linux
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows
"""

NAME = "miubrt"
MAINTAINER = "Kai Mühlbauer"
MAINTAINER_EMAIL = "kai.muehlbauer@uni-bonn.de"
DESCRIPTION = DOCLINES[0]
LONG_DESCRIPTION = "\n".join(DOCLINES[2:])
URL = "https://github.com/meteo-ubonn/miubrt"
DOWNLOAD_URL = "https://github.com/meteo-ubonn/miubrt"
LICENSE = "MIT"
CLASSIFIERS = list(filter(None, CLASSIFIERS.split("\n")))
PLATFORMS = ["Linux", "Mac OS-X", "Unix", "Windows"]

def setup_package():

    from setuptools import find_packages, setup

    with open("requirements.txt", "r") as f:
        INSTALL_REQUIRES = [rq for rq in f.read().split("\n") if rq != ""]

    # with open("requirements_optional.txt", "r") as f:
    #     OPTIONAL_REQUIRES = [rq for rq in f.read().split("\n") if rq != ""]

    with open("requirements_devel.txt", "r") as f:
        DEVEL_REQUIRES = [rq for rq in f.read().split("\n") if rq != ""]

    INSTALL_REQUIRES += OPTIONAL_REQUIRES

    metadata = dict(
        name=NAME,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        url=URL,
        download_url=DOWNLOAD_URL,
        license=LICENSE,
        classifiers=CLASSIFIERS,
        platforms=PLATFORMS,
        setup_requires=["setuptools_scm"],
        install_requires=INSTALL_REQUIRES,
        extras_require=dict(
            dev=DEVEL_REQUIRES,
        ),
        packages=find_packages(),
    )

    setup(**metadata)


if __name__ == "__main__":
    setup_package()
