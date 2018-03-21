#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""

"""
import os
import unittest

import numpy as np
from pygssi.lib import gssi_filelib
import hashlib


class TestRead(unittest.TestCase):

    def test_read_dzt(self):
        dst = gssi_filelib.read_dzt(os.path.split(__file__)[0] + '/data/test_dat1.DZT')
        self.assertTrue(dst.samp.shape == (4096, 500))
        self.assertTrue(dst.header.range == 1200)

    def test_read(self):
        gssi_filelib.read([os.path.split(__file__)[0] + '/data/test_dat1.DZT'], revs=[False])

    def test_cache(self):
        self.hashval = 'radar_' + hashlib.sha256(''.join([os.path.split(__file__)[0] + '/data/test_dat1.DZT']).encode('UTF-8')).hexdigest()
        stuff = gssi_filelib.read([os.path.split(__file__)[0] + '/data/test_dat1.DZT'], revs=[False], cache=True)
        stuff2 = gssi_filelib.read([os.path.split(__file__)[0] + '/data/test_dat1.DZT'], revs=[False])
        self.assertTrue(np.all(stuff[3][0].samp == stuff2[3][0].samp))
        os.remove(self.hashval)


class TestWrite(unittest.TestCase):

    def setUp(self):
        self.dzt = gssi_filelib.read_dzt(os.path.split(__file__)[0] + '/data/test_dat1.DZT')

    def testwrite(self):
        gssi_filelib.write_dzt(os.path.split(__file__)[0] + '/data/test_dat3.DZT', self.dzt)
        dz2 = gssi_filelib.read_dzt(os.path.split(__file__)[0] + '/data/test_dat3.DZT')
        self.assertTrue(np.all(self.dzt.samp == dz2.samp))

    def tearDown(self):
        os.remove(os.path.split(__file__)[0] + '/data/test_dat3.DZT')


if __name__ == '__main__':
    unittest.main()
