#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013)
#
# This file is part of the hveto python package.
#
# hveto is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hveto is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with hveto.  If not, see <http://www.gnu.org/licenses/>.

"""Setup the kp84 package
"""

from __future__ import print_function

import sys
if sys.version < '2.6':
    raise ImportError("Python versions older than 2.6 are not supported.")

import glob
import os.path

from setuptools import (setup, find_packages)

# set basic metadata
PACKAGENAME = 'kp84'
DISTNAME = 'kp84'
AUTHOR = 'Michael Coughlin'
AUTHOR_EMAIL = 'michael.coughlin@ligo.org'
LICENSE = 'GPLv3'

cmdclass = {}

# -- versioning ---------------------------------------------------------------

import versioneer
__version__ = versioneer.get_version()
cmdclass.update(versioneer.get_cmdclass())

# -- documentation ------------------------------------------------------------

# import sphinx commands
try:
    from sphinx.setup_command import BuildDoc
except ImportError:
    pass
else:
    cmdclass['build_sphinx'] = BuildDoc

# -- dependencies -------------------------------------------------------------

setup_requires = [
    'setuptools',
    'pytest-runner',
]
install_requires = [
    'astroml',
    'aplpy',
    'arrow',
    'astropy',
    'astroquery',
    'astroplan',
    'click',
    'Cython',
    'flask_caching',
    'flask_humanize',
    'flask_login',
    'flask_sqlalchemy',
    'flask_wtf',
    'extension_helpers',
    'h5py',
    'idna==2.8',
    'image_registration',
    'matplotlib',
    'numpy',
    'pandas',
    'passlib',
    'photutils',
    'pytest',
    'redis',
    'scipy',
    'scikit-image',
    'simplejson',
    'tables',
    'werkzeug',
    'WTForms',
    'WTForms-Alchemy',
    'WTForms-Components'
]
tests_require = [
    'pytest'
]
if sys.version_info < (2, 7):
    tests_require.append('unittest2')
extras_require = {
    'doc': [
        'sphinx',
        'numpydoc',
        'sphinx_rtd_theme',
        'sphinxcontrib_programoutput',
        'sphinxcontrib_epydoc',
    ],
}

# -- run setup ----------------------------------------------------------------

packagenames = find_packages()
scripts = glob.glob(os.path.join('bin', '*'))

setup(name=DISTNAME,
      provides=[PACKAGENAME],
      version=__version__,
      description=None,
      long_description=None,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      packages=packagenames,
      include_package_data=True,
      cmdclass=cmdclass,
      scripts=scripts,
      setup_requires=setup_requires,
      install_requires=install_requires,
      tests_require=tests_require,
      extras_require=extras_require,
      use_2to3=True,
      classifiers=[
          'Programming Language :: Python',
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Science/Research',
          'Natural Language :: English',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Physics',
          'Operating System :: POSIX',
          'Operating System :: Unix',
          'Operating System :: MacOS',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      ],
)
