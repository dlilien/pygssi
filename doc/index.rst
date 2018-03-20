.. pygssi documentation master file, created by
   sphinx-quickstart on Tue Oct  3 14:07:53 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pygssi's documentation!
==================================

`pygssi <http://www.github.com/dlilien/pygssi>`_ is an open-source library for opening GSSI radar files. This should give you access to all the information in the header (except the gain, which I have been unable to parse) and access to the raw returns. It can do some basic things like plot up those returns, and there is a (poor) tracing application. This library can also parse gga strings from a GPS attached to the SIR4000.

Requirements
------------

Python 2 or 3,
`numpy <http://www.scipy.org>`_,
`scipy <http://numpy.org>`_,
`pygeotools <http://www.github.com/dshean/pygeotools>`_.

Install
-------

You will first need to clone the repository using something like ``git clone http://www.github.com/dlilien/pygssi.git``. After going into the pygssi directory (``cd pygssi``), you should only need to run
``python setup.py install``
to gain access to the libraries and to get the executables on your path.

Testing
-------

Included are few basic tests. You will need to install `nose <http://nose.readthedocs.io/en/latest/>`_ to use them. They should be run with ``python setup.py test``

Documentation
=============

The main utilities to be aware of are:

`gssi_filelib <lib/gssi_filelib.html>`_, which contains the functions that you would want to use to import data for use in your own code

`gssilib <lib/gssilib.html>`_, which has some functions for filtering and plotting and overall processing

and `read_gssi.py <bin/read_gssi.html>`_, which will make a couple basic plots for you if all you need is a quick look. Further documentation can be found through the contents below.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   lib/index.rst
   bin/index.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
