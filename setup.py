from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext# as build_pyx
import numpy as np
from socket import gethostname
if gethostname() == 'LOULI':
        print gethostname()
        include_dirs = ['/home/ac/rsp/pyticles/',np.get_include(),'.']
        print include_dirs
else:
        include_dirs = ['/Users/acharles/masters/active/pyticles/trunk/',np.get_include(),'.']

setup(name = 'pyticles',
      packages=['pyticles'],
      package_dir={'pyticles':'.'},
      ext_modules=[
         Extension('pairsep', ['pairsep.pyx'],include_dirs=include_dirs),
         Extension('c_test', ['c_test.pyx'],include_dirs=include_dirs),
         Extension('c_forces', ['c_forces.pyx'],include_dirs=include_dirs),
         Extension('c_properties', ['c_properties.pyx'],include_dirs=include_dirs)
         ],
      cmdclass = { 'build_ext': build_ext })


