#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

# Define package info
name = "pycoMeth"
version = "0.2.9"
description = 'Python package for nanopore DNA methylation analysis downstream to Nanopolish'
with open("README.md", "r") as fh:
    long_description = fh.read()

# Collect info in a dictionnary for setup.py
setup(
    name = name,
    description = description,
    version = version,
    long_description = long_description,
    long_description_content_type="text/markdown",
    url = "https://github.com/a-slide/pycoMeth",
    author = 'Adrien Leger',
    author_email = 'aleg@ebi.ac.uk',
    license= "MIT",
    python_requires ='>=3.6',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3'],
    install_requires = [
        'numpy>=1.14.0',
        'tqdm>=4.23.4',
        "pandas>=0.25.1",
        "statsmodels>=0.10.1",
        "scipy>=1.3.1",
        "pyfaidx>=0.5.5.2"],
    packages = [name],
    entry_points = {'console_scripts': ['pycoMeth=pycoMeth.__main__:main']})
