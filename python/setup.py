from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from numpy import get_include
from glob import glob

extensions = [Extension (ext[:-4], [ext]) for ext in glob ('*.pyx')]

setup(
  name = 'vlasov-sheet',
  ext_modules = cythonize (extensions),
  include_dirs = [get_include ()],
)
