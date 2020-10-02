#!/usr/bin/env python
"""
setup.py file for reverse chirp example pycbc waveform plugin package
"""

from setuptools import Extension, setup, Command
from setuptools import find_packages

VERSION = '0.0.dev0'

setup (
    name = 'pycbc-waveform-SEOBNRE',
    version = VERSION,
    description = 'A SEOBNRE waveform plugin for PyCBC',
    long_description = open('descr.rst').read(),
    author = 'Yifan Wang',
    author_email = 'yifan.wang@aei.mpg.de',
    url = 'http://www.pycbc.org/',
    download_url = 'https://github.com/yi-fan-wang/pycbc-waveform-SEOBNRE/tarball/v%s' % VERSION,
    keywords = ['pycbc', 'waveform','signal processing', 'gravitational waves'],
    install_requires = ['pycbc'],
    py_modules = ['SEOBNRE'],
    entry_points = {"pycbc.waveform.td":"SEOBNRE = SEOBNRE:SEOBNRE_td",
                    "pycbc.waveform.fd":"SEOBNRE = SEOBNRE:SEOBNRE_fd"},
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
