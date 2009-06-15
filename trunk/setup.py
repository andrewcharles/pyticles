from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext as build_pyx

include_dirs = ['/Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages/numpy/core/include/']


setup(name = 'pyticles',
      packages=['pyticles'],
      package_dir={'pyticles':'.'},
      ext_modules=[
#         Extension('_integrator', ['integrator.pyx'],include_dirs=include_dirs),
#         Extension('_particles', ['particles_module.pyx'],include_dirs=include_dirs)
         Extension('c_test', ['c_test.pyx'],include_dirs=include_dirs),
         Extension('c_forces', ['c_forces.pyx'],include_dirs=include_dirs),
         Extension('c_properties', ['c_properties.pyx'],include_dirs=include_dirs)
         ],
      cmdclass = { 'build_ext': build_pyx })


