#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

# Define package info
name = "pycoMeth"
version = "0.4.13"
description = "DNA methylation analysis downstream to Nanopolish for Oxford Nanopore DNA sequencing datasets"
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
    license= "GPL",
    python_requires ='>=3.6',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3'],
    install_requires = [
        "numpy>=1.14.0",
        "scipy==1.4.1",
        "statsmodels>=0.11.1",
        "pandas>=1.0.3",
        "Jinja2>=2.11.1",
        "plotly==4.7.1",
        "pyfaidx>=0.5.8",
        "tqdm>=4.45.0",
        "colorlog>=4.1.0",
        "nbformat>=4.2.0",
        "kaleido"],
    packages = [name],
    package_dir = {name: name},
    package_data = {name: ['templates/*']},
    entry_points = {'console_scripts': ['pycoMeth=pycoMeth.__main__:main']})
