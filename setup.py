#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2016 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the MIT license.

"""
Created for compilation of fortran code
"""
import setuptools
from numpy.distutils.core import setup
setuptools

# I wrote a non-backwards compatible grid reader, so we need to check version
# Need to support python2 still because ryder does not have python3 c headers
if __name__ == '__main__':
    console_scripts = ['read_dzt.py=pygssi.bin.read_dzt:main', 'read_dzg.py=pygssi.bin.read_dzg:main', 'read_gssi.py=pygssi.bin.read_gssi:main', 'trace_gssi.py=pygssi.bin.trace_gssi:main']
    scripts = []
    setup(name='pygssi',
          version='0.1a1',
          description='Scripts for processing gssi radar data',
          url='http://github.com/dlilien/pygssi',
          author='David Lilien',
          author_email='dal22@uw.edu',
          license='MIT',
          scripts=scripts,
          console_scripts=console_scripts,
          entry_points={'console_scripts': console_scripts},
          install_requires=['numpy', 'scipy', 'matplotlib', 'pygeotools'],
          packages=['pygssi', 'pygssi.bin', 'pygssi.lib'],
          include_package_data=True,
          test_suite='nose.collector')
