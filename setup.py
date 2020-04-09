#!/usr/bin/env/python
"""Installation script
"""

import os
import numpy
from setuptools import setup, Extension
from setuptools.command.build_py import build_py as _build_py

LONG_DESCRIPTION = """ """

if os.environ.get("READTHEDOCS", False) == "True":
    INSTALL_REQUIRES = []
else:
    INSTALL_REQUIRES = ['numpy>=1.16.2']

if os.environ.get('PYTHON_INCLUDE_DIRS') is None:
    PYTHON_INCLUDE_DIRS = []
else:
    PYTHON_INCLUDE_DIRS = os.environ.get('PYTHON_INCLUDE_DIRS').split(':')


class build_py(_build_py):   
    def run(self):
        self.run_command("build_ext")
        return super().run()

sphereGIS = Extension(name='_sphereGIS', 
                    sources=['sphereGIS.i', 'sphereGIS.cpp'], 
                    depends=['sphereGIS.h'],
                    swig_opts=['-modern', '-c++'],
                    extra_compile_args=['-std=c++11'],
                    #libraries=['STARE'],
                    #library_dirs=STARE_LIB_DIRS, 
                    #include_dirs=INCLUDE_DIRS,   
                    language='c++')


# get all data dirs in the datasets module
data_files = []

setup(
    name='sphereGIS',
    version='0.0.1',
    description="",
    cmdclass={'build_py': build_py},
    long_description=LONG_DESCRIPTION,         
    py_modules = ['sphereGIS'],
    ext_modules=[sphereGIS],    
    python_requires=">=3.5",
    test_suite='tests',
    install_requires=INSTALL_REQUIRES
) 



