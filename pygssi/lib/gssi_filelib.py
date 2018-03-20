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

import os.path
import hashlib
import struct
import pickle
import numpy as np

from .gpslib import nmea_all_info
from .conversionlib import to_date
from .gpslib import kinematic_info
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

    # I think this implicitly detects if you are using a SIR3000 or SIR4000
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
    try:
        rh.Modify = to_date(modify_full)
    except ValueError:
        # I think the modification times are not actually good?
        rh.Modify = rh.Create

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
    data: :class:`~pygssi.lib.gpslib.NMEA`
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


def read(fns, revs, pickle_fn=None, elev_fn=None, t_srs='sps'):
    """Read in files (and try to get their GPS info too)

    Parameters
    ----------
    fns: iterable
        The files to read
    revs: iterable
        Reverse each file? (must be same length as fns)
    pickle_fn: str, optional
        Try to load pre-packaged results from this filename
    elev_fn: str, optional
        Load kinematic GPS data from this matlab file so that you don't use the crappy ones from the GPS on the radar
    t_srs: str, optional
        Load things into this projection (can convert to sps/EPSG:3031 or nps/EPSG:3413 from wgs84

    Returns
    -------
    gps_data: list of :class:`~pygssi.lib.gpslib.NMEA`
        GPS data for each of the input files
    """

    if len(fns) != len(revs):
        raise ValueError('Each file must have a bool to say whether to reverse it')
    hashval = 'radar_' + hashlib.sha256(''.join(fns).encode('UTF-8')).hexdigest()
    if pickle_fn is not None or os.path.exists(hashval):
        print('Loading pickled data')
        if pickle_fn is not None:
            hashval = pickle_fn
        (gps_data,
         stacked_data,
         total_lat,
         total_lon,
         total_dist,
         dzts,
         elev_list) = pickle.load(open(hashval, 'rb'))
    else:
        print('Loading data from DZT and DZG files')
        gps_data = [get_dzg_data(os.path.splitext(fn)[0] + '.DZG', t_srs, rev=rev) for fn, rev in zip(fns, revs)]

        # Now find the x coordinates for plotting
        for gpd in gps_data:
            gpd.set_proj('sps')
            gpd.get_dist()
        print('GPS data read')

        if elev_fn is not None:
            kinematic_data = kinematic_info(elev_fn)
            elev_list = [kinematic_data.query(gpd.x, gpd.y) for gpd in gps_data]
            print('Read in elevations from kinematic GPS')
        else:
            elev_list = None

        total_dist = np.hstack([gps_data[0].dist.flatten()[0:]] + [gps_data[i].dist.flatten()[1:] + np.sum([gps_data[j].dist.flatten()[-1] for j in range(i)]) for i in range(1, len(gps_data))])
        total_lat = np.hstack([gps_data[0].lat.flatten()[0:]] + [gps_data[i].lat.flatten()[1:] for i in range(1, len(gps_data))])
        total_lon = np.hstack([gps_data[0].lon.flatten()[0:]] + [gps_data[i].lon.flatten()[1:] for i in range(1, len(gps_data))])
        dzts = [read_dzt(fn, rev=rev) for fn, rev in zip(fns, revs)]
        check_headers(dzts)
        gps_stack_number = gps_data[0].scans[1] - gps_data[0].scans[0]

        # we are doing this next bit in two steps because of cutoff effects where the length of the gps and the stacked data don't match
        stack_data_list = [np.array([np.nanmean(dzts[j].samp[:, i * gps_stack_number:(i + 1) * gps_stack_number], axis=1) for i in range(dzts[j].samp.shape[1] // gps_stack_number)]).transpose() for j in range(len(dzts))]
        for dzt in dzts:
            dzt.samp = None
        for i in range(1, len(stack_data_list)):
            if stack_data_list[i].shape[1] == gps_data[i].dist.shape[0]:
                stack_data_list[i] = stack_data_list[i][:, :-1]
        stacked_data = np.hstack(stack_data_list)
        pickle.dump((gps_data, stacked_data, total_lat, total_lon, total_dist, dzts, elev_list), open(hashval, 'wb'))
    lldist = np.vstack((total_lon, total_lat, total_dist)).transpose()
    return gps_data, stacked_data, lldist, dzts, elev_list


class GssiError(Exception):
    pass
