from setuptools import setup, Extension

import os

cython_src= "src/bialignment"

USE_CYTHON = os.path.isfile(cython_src+".pyx")

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [Extension("bialignment", [cython_src+ext])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(cython_src+".pyx",
        annotate=False,
        compiler_directives = {'boundscheck': False,
            'language_level': 3})

with open("src/bialignment_nonpyx.py") as f:
    for line in f.readlines():
        if line.startswith('__version__'):
            VERSION = line.strip().split()[-1][1:-1]

with open("README.md", "r") as fh:
    long_description = fh.read()

CLASSIFIERS = """\
Development Status :: 5 - Production/Stable
Intended Audience :: Science/Research
License :: OSI Approved :: GNU General Public License v3 (GPLv3)
Programming Language :: Python :: 3
Programming Language :: Python :: 3 :: Only
Topic :: Scientific/Engineering
"""

setup(
    name = "bialign",
    version=VERSION,
    author='Sebastian Will',
    author_email='sebastian.will@polytechnique.edu',
    maintainer='Sebastian Will',
    maintainer_email='sebastian.will@polytechnique.edu',
    description='Bialignment of RNAs and proteins',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/s-will/BiAlign',
    classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
    package_dir={'':'src'},
    py_modules=['bialignment_nonpyx'],
    scripts=['src/bialign.py'],
    packages=[],
    zip_safe=False,
    ext_modules = extensions
)
