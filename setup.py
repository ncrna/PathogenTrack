#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import re
from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="PathogenTrack",
    version="0.2.0",
    author="Wei Zhang",
    author_email="admin@ncrna.net",
    include_package_data=True,
    description="A pipeline to identify pathogenic microorganisms from scRNA-seq raw data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ncrna/PathogenTrack",
    project_urls={
        "Quick Start": "https://github.com/ncrna/PathogenTrack/blob/master/doc/QUICK_START.md",
        "Bug Tracker": "https://github.com/ncrna/PathogenTrack/issues",
    },
    install_requires=[
        "biopython>=1.78",
        "python>=3.6",
    ],
    classifiers=[
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3',
    ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    entry_points={
        'console_scripts': [
            'PathogenTrack=PathogenTrack.PathogenTrack:main'
        ]
    },
    platforms=["all"],
    license='MIT',
)
