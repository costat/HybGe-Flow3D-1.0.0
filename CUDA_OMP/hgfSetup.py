#!/usr/bin/env python
"""
hgf.py for SWIG C++ libraries
"""
from distutils.core import *
from distutils      import sysconfig

import numpy
try:
  numpy_include = numpy.get_include()
except AttributeError:
  numpy_include = numpy.get_numpy_include()

hgf_module = Extension( '_hgf',
                       sources=['hgf.i', 'hgf.cpp', 'hgfMesh.cpp', 'hgfStokes.cpp', 'hgfPoreNetwork.cpp', 'hgfArrays.cpp', 'hgfBC.cpp', 'hgfIB.cpp', 'hgfPP.cpp'],
		       include_dirs = [numpy_include],
                       extra_compile_args=['-O3','-fopenmp'],
                       extra_link_args=['-lparalution'],
		       swig_opts=['-c++']
		      )

setup ( name = 'HybGe-Flow3d',
        version = '0.0.1',
        author = 'Timothy B. Costa',
	description = '2D/3D CFD software with Immersed Boundary',
	ext_modules = [hgf_module],
	py_modules = ["hgf"]
      )
