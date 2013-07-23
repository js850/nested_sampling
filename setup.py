import os
import numpy as np
from distutils.core import setup
from distutils.extension import Extension

## Numpy header files 
numpy_lib = os.path.split(np.__file__)[0] 
numpy_include = os.path.join(numpy_lib, 'core/include') 

setup(
    name="nested_sampling",
    version='0.1', 
    description="Python implementation of the nested sampling algorithm",
    url="https://github.com/js850/nested_sampling",
    #cmdclass = {'build_ext': build_ext},
    packages=["nested_sampling",
              "nested_sampling.src",
              "nested_sampling.models",
              "nested_sampling.utils",
             ],
    ext_modules= 
        [
            Extension("nested_sampling.src.cv_trapezoidal", ["nested_sampling/src/cv_trapezoidal.c", "nested_sampling/src/cv.c"],
                      include_dirs=[numpy_include],
                      extra_compile_args = ['-Wextra','-pedantic','-funroll-loops','-O3',],
                        ),
        ]
      )
