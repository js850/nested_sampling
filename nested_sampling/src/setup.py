import os
import numpy as np
from distutils.core import setup
from distutils.extension import Extension

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
            Extension("cv_trapezoidal", ["cv_trapezoidal.c", "cv.c"],
                      include_dirs=gslpy_include,
                      extra_compile_args = ['-Wextra','-pedantic','-funroll-loops','-O3',],
                      libraries=['gsl', 'gslcblas', 'm']
                        ),
        ]
      )
