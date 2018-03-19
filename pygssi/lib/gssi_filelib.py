#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dal22@uw.edu>
#
# Distributed under terms of the MIT license.

"""
Just define some classes and helper reader functions for gssi data
"""

import struct
import numpy as np
from .gpslib import nmea_all_info
from .conversionlib import to_date
from pygeotools.lib import geolib


class RH:
    """You probably need to read the GSSI docs on the c struct to understand this stuff"""
    tag = None
    data = None
    nsamp = None
    bits = None
    bytes = None
    us_dattype = None
    s_dattype = None
    rgain = None
    nrgain = None
    checksum = None
    antname = None

    def __str__(self):
        return 'rgain: {:d}, nrgain {:d}'.format(self.rgain, self.nrgain)

    def __repr__(self):
        return self.__str__()


class DZT:
    """
    A DZT file

    Attributes
    ----------
    header:  :class:`~pygssi.lib.gssi_filelib.RH`
        The header of the file, with lot of useful (maybe) info
    samp: numpy.ndarray
        An array containing the actual returns (not normalized, so won't plot well as-is
    """
    header = None
    samp = None

    def __init__(self, header, sample):

        self.header = header
        self.samp = sample


def read_dzt(fn, rev=False):
    """
    Read a dzt file

    Parameters
    ----------
    fn: str
        File to parse
    rev: bool, optional
        Reverse the input array (useful for concatenating files)

    Returns
    -------
    data: :class:`~pygssi.lib.gssi_filelib.DZT`
    """

    rh = RH()
    with open(fn, 'rb') as fid:
        lines = fid.read()
    rh.tag = struct.unpack('<H', lines[0:2])[0]
    rh.data = struct.unpack('<H', lines[2:4])[0]
    rh.nsamp = struct.unpack('<H', lines[4:6])[0]
    rh.bits = struct.unpack('<H', lines[6:8])[0]
    rh.bytes = rh.bits // 8

    # I think this implicitly detects if you are using a SIR2000 or SIR3000
    if rh.bits == 32:
        rh.us_dattype = 'I'
    elif rh.bits == 16:
        rh.us_dattype = 'H'
    if rh.bits == 32:
        rh.s_dattype = 'i'
    elif rh.bits == 16:
        rh.s_dattype = 'h'

    rh.zero = struct.unpack('<h', lines[8:10])[0]
    rh.sps = struct.unpack('<f', lines[10:14])[0]
    rh.spm = struct.unpack('<f', lines[14:18])[0]
    rh.mpm = struct.unpack('<f', lines[18:22])[0]
    rh.position = struct.unpack('<f', lines[22:26])[0]
    rh.range = struct.unpack('<f', lines[26:30])[0]

    rh.npass = struct.unpack('<h', lines[30:32])[0]

    create_full = struct.unpack('<4s', lines[32:36])[0]
    modify_full = struct.unpack('<4s', lines[36:40])[0]

    rh.Create = to_date(create_full)
    rh.Modify = to_date(modify_full)

    rh.rgain = struct.unpack('<H', lines[40:42])[0]
    rh.nrgain = struct.unpack('<H', lines[42:44])[0] + 2
    rh.text = struct.unpack('<H', lines[44:46])[0]
    rh.ntext = struct.unpack('<H', lines[46:48])[0]
    rh.proc = struct.unpack('<H', lines[48:50])[0]
    rh.nproc = struct.unpack('<H', lines[50:52])[0]
    rh.nchan = struct.unpack('<H', lines[52:54])[0]

    rh.epsr = struct.unpack('<f', lines[54:58])[0]
    rh.top = struct.unpack('<f', lines[58:62])[0]
    rh.depth = struct.unpack('<f', lines[62:66])[0]

    rh.reserved = struct.unpack('<31c', lines[66:97])
    rh.dtype = struct.unpack('<c', lines[97:98])[0]
    rh.antname = struct.unpack('<14c', lines[98:112])

    rh.chanmask = struct.unpack('<H', lines[112:114])[0]
    rh.name = struct.unpack('<12c', lines[114:126])
    rh.chksum = struct.unpack('<H', lines[126:128])[0]

    rh.breaks = struct.unpack('<H', lines[rh.rgain:rh.rgain + 2])[0]
    rh.Gainpoints = np.array(struct.unpack('<{:d}i'.format(rh.nrgain), lines[rh.rgain + 2:rh.rgain + 2 + 4 * (rh.nrgain)]))
    rh.Gain = 0
    if rh.ntext != 0:
        rh.comments = struct.unpack('<{:d}s'.format(rh.ntext), lines[130 + 2 * rh.Gain: 130 + rh.bytes * rh.Gain + rh.ntext])[0]
    else:
        rh.comments = ''
    if rh.nproc != 0:
        rh.proccessing = struct.unpack('<{:d}s'.format(rh.nproc), lines[130 + rh.bytes * rh.Gain + rh.ntext:130 + rh.bytes * rh.Gain + rh.ntext + rh.nproc])[0]
    else:
        rh.proc = ''

    # d = np.array(struct.unpack('<{:d}'.format((len(lines) - 1024) // 2) + dattype, lines[1024:]))
    d = np.array(struct.unpack('<{:d}'.format((len(lines) - 36 * 4096) // rh.bytes) + rh.us_dattype, lines[36 * 4096:]))
    d = d.reshape((rh.nsamp, -1), order='F')
    d[0, :] = d[2, :]
    d[1, :] = d[2, :]
    d = d + rh.zero
    if rev:
        d = np.fliplr(d)
    dat = DZT(rh, d)
    return dat


def get_dzg_data(fn, t_srs='sps', rev=False):
    """Read GPS data associated with a GSSI sir4000 file.

    Parameters
    ----------
    fn: str
        A dzg file with ggis and gga strings.
    t_srs: str, optional
        Target coordinate reference. Default lat/lon (wgs84)
    rev: bool, optional
        Reverse the points in this file (used for concatenating radar files). Default False.
    
    Returns
    -------
    data: :class:`~pygssi.lib.gpslib.nmea_info`
    """

    with open(fn) as f:
        lines = f.readlines()
    ggis = lines[::3]
    gga = lines[1::3]
    data = nmea_all_info(gga)
    data.scans = np.array(list(map(lambda x: int(x.split(',')[1]), ggis)))
    if rev:
        data.rev()
    data.get_all()

    if t_srs in ['sps', 'EPSG:3031']:
        data.x, data.y, _ = geolib.ll2sps(data.lon, data.lat, data.z)
    elif t_srs in ['nps', 'EPSG:3413']:
        data.x, data.y, _ = geolib.ll2nps(data.lon, data.lat, data.z)
    else:
        print('t_srs not recognized, defaulting to lat/lon (wgs84)')
        data.x, data.y = data.lon, data.lat
    return data


def check_headers(dzts):
    """This function should check that the headers have the same number of channels, bytes, etc. and raise an exception if not"""
    if False:
        raise GssiError


class GssiError(Exception):
    pass
