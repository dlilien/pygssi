#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""

"""
import os.path
import unittest

from pygssi.lib import gssilib


class TestProcessRadar(unittest.TestCase):

    def test_bare(self):
        gssilib.process_radar(os.path.split(__file__)[0] + '/data/test_dat1.DZT', cache=False)
        self.assertTrue(True)

    def test_rev(self):
        gssilib.process_radar(os.path.split(__file__)[0] + '/data/test_dat1.DZT', cache=False, rev_list=True)
        self.assertTrue(True)

    def test_two(self):
        gssilib.process_radar([os.path.split(__file__)[0] + '/data/test_dat1.DZT', os.path.split(__file__)[0] + '/data/test_dat2.DZT'], cache=False, rev_list=True)
        self.assertTrue(True)

    def test_tworev(self):
        gssilib.process_radar([os.path.split(__file__)[0] + '/data/test_dat1.DZT', os.path.split(__file__)[0] + '/data/test_dat2.DZT'], cache=False, rev_list=[False, True])
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
