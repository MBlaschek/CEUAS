from setuptools import setup

setup(name='CEUAS',
      version='0.1',
      description='Copernicus Climate Change Early Upper Air Service',
      url='https://github.com/MBlaschek/CEUAS',
      author='MB',
      author_email='michael.blaschek@univie.ac.at',
      license='UNIVIE GNU GPL',
      packages=['CEUAS', 'doc'],
      install_requires=['numpy', 'pandas', 'xarray', 'numba', 'requests'],
      zip_safe=False)
