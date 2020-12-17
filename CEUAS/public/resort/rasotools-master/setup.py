import re

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt', 'r') as f:
    required = [re.sub(r'==.*', '', i) for i in f.read().splitlines()]

# date based versioning [YEAR].[MONTH]
setuptools.setup(
    name='rasotools',
    version='19.12',
    description='RASO tools',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/MBlaschek/rasotools',
    author='MB',
    author_email='michael.blaschek@univie.ac.at',
    license='MIT',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
    ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=required,
    scripts=['bin/rsncnfo'],
    entry_points={
        'console_scripts': [
            'rscf = rasotools.fun.cfconvention:main'
        ]
    },
    python_requires='>=3'
)

# Build and upload
# python3 setup.py bdist_wheel
# twine upload dist/* -r pypitest --verbose
# make requirements.txt
# pipreqs rasotools  -> rasotools/requirements.txt
