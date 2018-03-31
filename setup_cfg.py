#! /usr/bin/env python

# System imports
from distutils.core import *
from distutils      import sysconfig

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the gsl include and library directories.
GSL_INCLUDE = '/opt/local/include'
GSL_LIB = '/opt/local/lib'

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    NUMPY_INCLUDE = numpy.get_include()
except AttributeError:
    NUMPY_INCLUDE = numpy.get_numpy_include()

source_files = ['bstrap_bc_wrap.c','characterization.c','prob.c','prob_util.c','nsb.c']
source_files = [('src/' + x) for x in source_files]

# range extension module
bstrap_bc = Extension('thoth._bstrap_bc',
                   sources = source_files,
                   library_dirs = [GSL_LIB],
                   libraries = ['gsl', 'gslcblas'],
                   include_dirs = [NUMPY_INCLUDE, GSL_INCLUDE],
		   extra_link_args = ['-fpic'],
		   runtime_library_dirs = [GSL_LIB])

# NumyTypemapTests setup
setup(  name        = "THOTH",
        description = "demo",
        author      = "Simon DeDeo",
        version     = "1.0",
        package_dir = {'thoth':'src'},
        packages = ['thoth'],
        ext_modules = [bstrap_bc],
        py_modules = ["src/thoth"]
        )
