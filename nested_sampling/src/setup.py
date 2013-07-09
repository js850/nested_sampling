import os
import numpy as np
from distutils.core import setup
from distutils.extension import Extension
#from Cython.Distutils import build_ext

## Python include 
#py_include = get_python_inc() 

gsl_include = '/usr/local/include' 

## Numpy header files 
numpy_lib = os.path.split(np.__file__)[0] 
numpy_include = os.path.join(numpy_lib, 'core/include') 

## Gslpy include 
gslpy_include = [gsl_include, numpy_include] 



setup(
    #cmdclass = {'build_ext': build_ext},
    ext_modules = 
        [
            Extension("minima_sampling", ["weighted_pick.c"]),
            Extension("runmc", ["runmc.c", "mc.c", "lj.c"],
                      include_dirs=gslpy_include,
                      libraries=['gsl', 'gslcblas', 'm']
                        ),
        ]
)
