pygssi
======
[![Build Status](https://travis-ci.org/dlilien/pygssi.svg?branch=master)](https://travis-ci.org/dlilien/pygssi)

An open-source library for opening GSSI radar files. This should give you access to all the information in the header (except the gain, which I have been unable to parse) and access to the raw returns. This library can also parse gga strings from a GPS attached to the SIR4000. [Read the docs.](https://dlilien.github.io/pygssi)

Requirements
------------

Python 2 or 3,
[numpy](http://www.scipy.org),
[scipy](http://numpy.org),
[pygeotools](http://www.github.com/dshean/pygeotools).

Install
-------

You should only need to run
``python setup.py install``
to gain access to the libraries and to get the executables on your path.

Documentation
-------------

A fairly complete description of the functionality of this package can be found [here](https://dlilien.github.io/pygssi).
