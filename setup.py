import os
from setuptools import setup, find_packages
execfile(os.path.dirname(os.path.realpath(__file__)) + "/version.py")
setup(name='pavfinder',
      version=__version__,
      author='Readman Chiu',
      packages=find_packages()
      )
