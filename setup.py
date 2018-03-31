#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""
# System imports
import os.path
import platform
from setuptools import setup, find_packages, Extension

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the gsl include and library directories.
GSL_INCLUDE = os.environ.get('GSL_INCLUDE')
GSL_LIB = os.environ.get('GSL_LIB')

# Set GSL_INCLUDE and GSL_LIB manually here, if necessary:
# GSL_INCLUDE = "include/"
# GSL_LIB = "lib/"


def get_gsl_paths(_gsl_include=None, _gsl_lib=None):
    system = platform.system()
    distro = None
    gsl_include, gsl_lib = None, None
    if system == "Darwin":
        # MacPorts: sudo port install gsl +universal
        gsl_include = '/opt/local/include'
        gsl_lib = '/opt/local/lib'
    elif system == "Linux":
        # distro = platform.linux_distribution()
            # Fedora: dnf install -y gsl gsl-devel
            # Ubuntu: apt
        gsl_include = '/usr/include/gsl'
        gsl_lib = '/usr/lib'
    else:
        # raise Exception()  # the default gcc paths may be fine
        print((
            "Unknown platform %r (%r):"
            "please define GSL_INCLUDE and GSL_LIB manually in setup.py") %
                (system, distro))
    return _gsl_include or gsl_include, _gsl_lib or gsl_lib


GSL_INCLUDE, GSL_LIB = get_gsl_paths(GSL_INCLUDE, GSL_LIB)

include_dirs = []
library_dirs = []
runtime_library_dirs = []

# Recognize PREFIX or CONDA_PREFIX ./include and ./lib
prefix = os.environ.get('PREFIX', os.environ.get('CONDA_PREFIX'))
if prefix is not None:
    include_dirs.append(os.path.join(prefix, 'include'))
    library_dirs.append(os.path.join(prefix, 'lib'))

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    NUMPY_INCLUDE = numpy.get_include()
except AttributeError:
    NUMPY_INCLUDE = numpy.get_numpy_include()

library_dirs.append(GSL_LIB) if GSL_LIB else None
include_dirs.append(NUMPY_INCLUDE)
include_dirs.append(GSL_INCLUDE) if GSL_INCLUDE else None
runtime_library_dirs.append(GSL_LIB) if GSL_LIB else None

source_files = [
    'bstrap_bc_wrap.c', 'characterization.c', 'prob.c', 'prob_util.c', 'nsb.c']
source_files = [os.path.join('thoth', x) for x in source_files]

# range extension module
bstrap_bc = Extension(
    'thoth._bstrap_bc',
    sources=source_files,
    library_dirs=library_dirs,
    libraries=['gsl', 'gslcblas'],
    include_dirs=include_dirs,
    extra_link_args=['-fpic'],
    runtime_library_dirs=runtime_library_dirs)


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['numpy']

setup_requirements = ['numpy']

test_requirements = []

setup(
    author="Simon DeDeo",
    author_email='sdedeo@andrew.cmu.edu',
    classifiers=[
        #'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        # 'Programming Language :: Python :: 3',
        # 'Programming Language :: Python :: 3.4',
        # 'Programming Language :: Python :: 3.5',
        # 'Programming Language :: Python :: 3.6',
    ],
    description=(
        "a python package for the efficient estimation of "
        "information-theoretic quantities from empirical data"),
    install_requires=requirements,
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='thoth entropy informationtheory probability',
    name='thoth',
    packages=find_packages(include=['thoth']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/simondedeo/thoth',
    version='1.0.0',
    zip_safe=False,
    ext_modules=[bstrap_bc],
)
