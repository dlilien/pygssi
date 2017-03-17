#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2016 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""
Some nmea conversion
"""
from pygeotools.lib import geolib
import numpy as np


class nmea_info:
    all_data = None
    lat = None
    lon = None
    qual = None
    sats = None
    x = None
    y = None
    z = None
    geo_offset = None
    times = None
    scans = None

    def rev(self):
        self.all_data = np.flipud(self.all_data)

    def get_all(self):
        self.glat()
        self.glon()
        self.gqual()
        self.gsats()
        self.gz()
        self.ggeo_offset()
        self.gtimes()

    def glat(self):
        self.lat = self.all_data[:, 2] * ((self.all_data[:, 1] - self.all_data[:, 1] % 100) / 100 + (self.all_data[:, 1] % 100) / 60)
        if self.y is None:
            self.y = self.lat
        return self.lat

    def glon(self):
        self.lon = self.all_data[:, 4] * ((self.all_data[:, 3] - self.all_data[:, 3] % 100) / 100 + (self.all_data[:, 3] % 100) / 60)
        if self.x is None:
            self.x = self.lon
        return self.lon

    def gqual(self):
        self.qual = self.all_data[:, 5]
        return self.qual

    def gsats(self):
        self.sats = self.all_data[:, 6]
        return self.sats

    def gz(self):
        self.z = self.all_data[:, 7]
        return self.z

    def ggeo_offset(self):
        self.geo_offset = self.all_data[:, 8]
        return self.geo_offset

    def gtimes(self):
        self.times = self.all_data[:, 9]
        return self.times

    def set_proj(self, proj_name='sps'):
        if self.lat is None:
            self.glat()
        if self.lon is None:
            self.glon()
        if self.z is None:
            self.gz()

        if proj_name in ['ll', 'wgs84']:
            return None
        elif proj_name in ['EPSG:3031', 'sps']:
            self.x, self.y, self.z = geolib.ll2sps(self.lon, self.lat, self.z)
            return self.x, self.y
        raise ValueError('Unrecognized coordinate system')

    def get_dist(self):
        if self.x is None:
            self.glat()
        if self.y is None:
            self.glon()
        self.dist = np.zeros((len(self.y), 1))
        for i in range(1, len(self.dist)):
            self.dist[i] = self.dist[i - 1] + np.sqrt((self.x[i] - self.x[i - 1]) ** 2.0 + (self.y[i] - self.y[i - 1]) ** 2.0)


def nmea_to_ll(list_of_sentences):
    def _gga_sentence_to_ll(sentence):
        vals = sentence.split(',')
        lat = float(vals[2])
        lon = float(vals[4])
        return (lon, lat)

    if list_of_sentences[0].split(',')[0] == '$GPGGA':
        return np.array([_gga_sentence_to_ll(sentence) for sentence in list_of_sentences]) / 100.0
    else:
        raise ValueError('I can only do gga sentences right now')


def nmea_all_info(list_of_sentences):
    """Return an object with the nmea info"""
    def _gga_sentence_split(sentence):
        all = sentence.split(',')
        numbers = list(map(float, all[1:3] + [1] + [all[4]] + [1] + all[6:10] + [all[11]]))
        if all[3] == 'S':
            numbers[2] = -1
        if all[5] == 'W':
            numbers[4] = -1
        return numbers

    if list_of_sentences[0].split(',')[0] == '$GPGGA':
        data = nmea_info()
        data.all_data = np.array([_gga_sentence_split(sentence) for sentence in list_of_sentences])
        return data
    else:
        print(list_of_sentences[0].split(',')[0])
        raise ValueError('I can only do gga sentences right now')
