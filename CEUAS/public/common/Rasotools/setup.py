from distutils.core import setup

setup(
    name='Rasotools',
    version='0.1',
    author='Leopold Haimberger',
    author_email='leopold.haimberger@univie.ac.at',
    packages=['rasotools', 'rasotools.test'],
    scripts=['bin/iscript.py'],
    url='http://pypi.python.org/pypi/Rasotools/',
    license='LICENSE.txt',
    description='Utilities for manipulating netcdf CF1.0 radiosonde data',
    long_description=open('README.txt').read(),
    install_requires=[
        "Numpy >= 1.8",
        "numexpr >= 2.4",
    ],
)