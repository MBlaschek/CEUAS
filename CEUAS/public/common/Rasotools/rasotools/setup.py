from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'anomaly calculations',
  ext_modules = cythonize("anomaly.pyx"),
)

