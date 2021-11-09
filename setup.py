from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("bialignment.pyx",
        annotate=False,
        compiler_directives = {'boundscheck': False,
            'language_level': 3})
    
)

