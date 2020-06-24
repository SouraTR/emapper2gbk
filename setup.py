#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    "biopython",
    "coloredlogs",
    "gffutils",
    "pandas",
    "pronto",
    "requests"]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]

description = "Build .gbk files starting from eggnog annotation files and genomes (fasta)"
# In case of long description
# description +=

setup(
    author="Clémence Frioux",
    author_email='clemence.frioux@inria.fr',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description=description,
    entry_points={
        'console_scripts': [
            'emapper2gbk=emapper2gbk.__main__:cli',
        ],
    },
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme + history,
    include_package_data=True,
    keywords='emapper2gbk',
    name='emapper2gbk',
    packages=find_packages(include=['emapper2gbk', 'emapper2gbk.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/cfrioux/emapper2gbk',
    zip_safe=False,
)
