#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2016 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""
Class definitions and helper functions for gssi info
"""
import struct
import numpy as np
from .gpslib import nmea_all_info, kinematic_info
import os
import codecs
import pickle
import hashlib
from math import atan2, degrees


import matplotlib.pyplot as plt
from cycler import cycler
from scipy import signal
from scipy.stats import linregress
from scipy.io import loadmat

from pygeotools.lib import geolib

# No gray
plt.rc('axes', prop_cycle=(cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#17becf']) + cycler('linestyle', ['solid'] * 9)))


def labelLine(line, x, label=None, align=True, usex=True, bg=True, **kwargs):
    """Label a line. Can do it based on x or y position. This is in here for labeling lat/lon

    Parameters
    ----------
    line: matplotlib.pyplot.Line2D
        The line you want to label
    x: float
        The x or y position to label (depending on usex)
    label: str, optional
        What the label should say
    usex: bool, optional
        Label based on x position. If false second argument gives the y position. Default True
    kwargs: passed on to matplotlib.pyplot.text
    """

    ax = line.axes
    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if ((x < xdata[0]) or (x > xdata[-1])) and usex:
        xdata = xdata[::-1]
        reverse_x = True
        if ((x < xdata[0]) or (x > xdata[-1])):
            # print('X {:f} outside bounds {:f} {:f}'.format(x, xdata[0], xdata[-1]))
            return
    else:
        reverse_x = False

    if ((x < ydata[0]) or (x > ydata[-1])) and not usex:
        ydata = ydata[::-1]
        reverse_y = True
        if ((x < ydata[0]) or (x > ydata[-1])):
            # print('Y {:f} outside bounds {:f} {:f}'.format(x, ydata[0], ydata[-1]))
            return
    else:
        reverse_y = False

    # Find corresponding y co-ordinate and angle of the
    if usex:
        if reverse_x:
            for i in range(len(xdata)):
                if x > xdata[i]:
                    break
            xdata = xdata[::-1]
        else:
            for i in range(len(xdata)):
                if x < xdata[i]:
                    break

        y = ydata[i - 1] + (ydata[i] - ydata[i - 1]) * \
            (x - xdata[i - 1]) / (xdata[i] - xdata[i - 1])
    else:
        y = x
        if reverse_y:
            for i in range(len(ydata)):
                if y > ydata[i]:
                    break
            ydata = ydata[::-1]
        else:
            for i in range(len(ydata)):
                if y < ydata[i]:
                    break
        x = xdata[i - 1] + (xdata[i] - xdata[i - 1]) * (y - ydata[i - 1]) / (ydata[i] - ydata[i - 1])
        if np.isnan(x):
            x = xdata[i]

    if not label:
        label = line.get_label()

    if align:
        # Compute the slope
        dx = xdata[i] - xdata[i - 1]
        dy = ydata[i] - ydata[i - 1]
        if dx < 0 and dy > 0:
            dx = -dx
            dy = -dy
        ang = degrees(atan2(dy, dx))
        if np.abs(90 - ang) < 1.0e-3:
            ang = -ang

        # Transform to screen co-ordinates
        pt = np.array([x, y]).reshape((1, 2))
        trans_angle = ax.transData.transform_angles(np.array((ang,)), pt)[0]

    else:
        trans_angle = 0

    # Set a bunch of keyword arguments
    if 'color' not in kwargs:
        kwargs['color'] = line.get_color()

    if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
        kwargs['ha'] = 'left'

    if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
        kwargs['va'] = 'center'

    if 'backgroundcolor' not in kwargs:
        kwargs['backgroundcolor'] = ax.get_facecolor()

    if 'clip_on' not in kwargs:
        kwargs['clip_on'] = True

    if 'zorder' not in kwargs:
        kwargs['zorder'] = 2.5

    return ax.text(x, y, label, rotation=trans_angle, **kwargs)


class DZT:
    header = None
    samp = None

    def __init__(self, header, sample):
        self.header = header
        self.samp = sample


class time:
    sec2 = None
    minute = None
    hour = None
    day = None
    month = None
    year = None


class RH:
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


def bits(bytes):
    for b in bytes:
        for i in range(8):
            yield (b >> i) & 1


def to_date(bin, le=True):
    a = time()
    bit = [b for b in bits(bin)]
    a.sec2 = bit_to_int(bit[0:5])
    a.minute = bit_to_int(bit[5:11])
    a.hour = bit_to_int(bit[11:16])
    a.day = bit_to_int(bit[16:21])
    a.month = bit_to_int(bit[21:25])
    a.year = bit_to_int(bit[25:32])
    return a


def bit_to_int(bits):
    return sum([(2 ** i) * bit for i, bit in enumerate(bits)])


def read_dzt(fn, rev=False):
    rh = RH()
    with open(fn, 'rb') as fid:
        lines = fid.read()
    rh.tag = struct.unpack('<H', lines[0:2])[0]
    rh.data = struct.unpack('<H', lines[2:4])[0]
    rh.nsamp = struct.unpack('<H', lines[4:6])[0]
    rh.bits = struct.unpack('<H', lines[6:8])[0]
    rh.bytes = rh.bits // 8
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
    return data


def check_headers(dzts):
    """This function should check that the headers have the same number of channels, bytes, etc. and raise an exception if not"""
    pass


def tt_to_m_variable(diel_array, tt):
    # I am not sure if there is a smart way to do this, so instead I'm going to figure out ho
    current_layer = 0
    current_dist = 0.
    remaining_inc = tt

    # this is complicated in case we get a sparsely sampled file where layers get skipped
    while current_layer < diel_array.shape[0]:
        if current_dist + remaining_inc * 1.0e-9 * 3.0e8 / 2. / np.sqrt(diel_array[current_layer, 1]) <= diel_array[current_layer, 0]:
            y = current_dist + remaining_inc * 1.0e-9 * 3.0e8 / 2. / np.sqrt(diel_array[current_layer, 1])
            break
        else:
            dist_this_layer = diel_array[current_layer, 0] - current_dist
            time_for_this_layer = dist_this_layer / 1.0e-9 / 3.0e8 * 2. * np.sqrt(diel_array[current_layer, 1])
            remaining_inc -= time_for_this_layer
            current_dist = diel_array[current_layer, 0]
            current_layer += 1
    return y


def tt_to_m_variable_arr(diel_array, tt_arr):
    y = np.zeros(tt_arr.shape)
    # I am not sure if there is a smart way to do this, so instead I'm going to figure out ho
    current_layer = 0
    current_dist = 0.

    # this is complicated in case we get a sparsely sampled file where layers get skipped
    while current_layer < diel_array.shape[0]:
        bool_arr = np.logical_and(y == 0, current_dist + tt_arr * 1.0e-9 * 3.0e8 / 2. / np.sqrt(diel_array[current_layer, 1]) <= diel_array[current_layer, 0])
        y[bool_arr] = (current_dist + tt_arr * 1.0e-9 * 3.0e8 / 2. / np.sqrt(diel_array[current_layer, 1]))[bool_arr]

        dist_this_layer = diel_array[current_layer, 0] - current_dist
        time_for_this_layer = dist_this_layer / 1.0e-9 / 3.0e8 * 2. * np.sqrt(diel_array[current_layer, 1])
        tt_arr -= time_for_this_layer
        current_dist = diel_array[current_layer, 0]
        current_layer += 1
    return y


def load_dielf(fn):
    with open(fn) as fin:
        lines = fin.readlines()

    # pad this with some values so we don't go out of bounds
    lines = list(map(lambda x: list(map(float, x.replace('\ufeff', '').rstrip('\n\r').split(','))), lines))
    diels = np.array([[0.0, lines[0][1]]] + lines + [[np.inf, lines[-1][1]]])
    return diels


def process_radar(fns,
                  rev_list,
                  t_srs='sps',
                  elev_fn=None,
                  pickle_fn=None,
                  axin=None,
                  gp=None,
                  dielf=None,
                  diel=None,
                  layers=None,
                  slope=False,
                  markersize=5,
                  label=False,
                  lw=1,
                  plotdata=True,
                  with_elevs=True,
                  xoff=0.0,
                  plot_layer=None,
                  bold_layer=None):
    hashval = hashlib.sha256(''.join(fns).encode('UTF-8')).hexdigest()
    if pickle_fn is not None or os.path.exists(hashval):
        print('Loading pickled data')
        if pickle_fn is not None:
            hashval = pickle_fn
        (gps_data,
         stacked_data,
         kinematic_data,
         total_lat,
         total_lon,
         total_dist,
         dzts,
         elev_list) = pickle.load(open(hashval, 'rb'))
    else:
        print('Loading data from DZT and DZG files')
        gps_data = [get_dzg_data(os.path.splitext(fn)[0] + '.DZG', t_srs, rev=rev) for fn, rev in zip(fns, rev_list)]
        # Now find the x coordinates for plotting
        for gpd in gps_data:
            gpd.set_proj('sps')
            gpd.get_dist()

        if elev_fn is not None:
            kinematic_data = kinematic_info(elev_fn)
            elev_list = [kinematic_data.query(gpd.x, gpd.y) for gpd in gps_data]
            print('Read in elevations')

        total_dist = np.hstack([gps_data[0].dist.flatten()[0:]] + [gps_data[i].dist.flatten()[1:] + np.sum([gps_data[j].dist.flatten()[-1] for j in range(i)]) for i in range(1, len(gps_data))])
        total_lat = np.hstack([gps_data[0].lat.flatten()[0:]] + [gps_data[i].lat.flatten()[1:] for i in range(1, len(gps_data))])
        total_lon = np.hstack([gps_data[0].lon.flatten()[0:]] + [gps_data[i].lon.flatten()[1:] for i in range(1, len(gps_data))])
        dzts = [read_dzt(fn, rev=rev) for fn, rev in zip(fns, rev_list)]
        # dzts = pickle.load(open('pickled_dzt', 'rb'))
        check_headers(dzts)
        # gps_data = pickle.load(open('pickled_gps', 'rb'))
        gps_stack_number = gps_data[0].scans[1] - gps_data[0].scans[0]

        # we are doing this next bit in two steps because of cutoff effects where the length of the gps and the stacked data don't match
        stack_data_list = [np.array([np.nanmean(dzts[j].samp[:, i * gps_stack_number:(i + 1) * gps_stack_number], axis=1) for i in range(dzts[j].samp.shape[1] // gps_stack_number)]).transpose() for j in range(len(dzts))]
        for dzt in dzts:
            dzt.samp = None
        for i in range(1, len(stack_data_list)):
            if stack_data_list[i].shape[1] == gps_data[i].dist.shape[0]:
                stack_data_list[i] = stack_data_list[i][:, :-1]
        stacked_data = np.hstack(stack_data_list)

        pickle.dump((gps_data, stacked_data, kinematic_data, total_lat, total_lon, total_dist, dzts, elev_list), open(hashval, 'wb'))

    total_dist += xoff

    if elev_fn is not None:
        elev = np.concatenate(elev_list).flatten()[:len(total_dist)]
    else:
        elev = np.zeros(total_dist.shape).flatten()

    stacked_data, zero_ind = zero_surface(stacked_data)

    if gp is not None:
        stacked_data = gained_decibels(gp, stacked_data)
    else:
        stacked_data = data_to_db(stacked_data)

    if dielf is not None:
        diels = load_dielf(dielf)
        
        layer_time_increment = dzts[0].header.range / dzts[0].header.nsamp

        # I am not sure if there is a smart way to do this, so instead I'm going to figure out ho
        y = np.zeros((stacked_data.shape[0], ))
        current_layer = 0
        current_dist = 0.
        y[0] = current_dist
        for i in range(1, len(y)):
            remaining_inc = layer_time_increment

            # this is complicated in case we get a sparsely sampled file where layers get skipped
            while current_layer < diels.shape[0]:
                if current_dist + remaining_inc * 1.0e-9 * 3.0e8 / 2. / np.sqrt(diels[current_layer, 1]) <= diels[current_layer, 0]:
                    y[i] = current_dist + remaining_inc * 1.0e-9 * 3.0e8 / 2. / np.sqrt(diels[current_layer, 1])
                    current_dist = y[i]
                    break
                else:
                    dist_this_layer = diels[current_layer, 0] - current_dist
                    time_for_this_layer = dist_this_layer / 1.0e-9 / 3.0e8 * 2. * np.sqrt(diels[current_layer, 1])
                    remaining_inc -= time_for_this_layer
                    current_dist = diels[current_layer, 0]
                    current_layer += 1
        y = np.flipud(y)
    elif diel is not None:
        # Find the travel time per pixel
        incremental_travel_time = dzts[0].header.range * 1.0e-9 * 3.0e8 / np.sqrt(diel) / 2. / dzts[0].header.nsamp
        y = np.flipud(np.arange(stacked_data.shape[0])) * incremental_travel_time
    else:
        y = np.flipud(np.arange(stacked_data.shape[0]))

    if axin is None:
        fig, ax = plot_bens_radar(stacked_data, x=total_dist / 1000.0, y=y, elev=elev)
        fig2, ax2 = plot_bens_radar(stacked_data, x=total_dist / 1000.0, y=y, elev=None)
    else:
        if plotdata:
            if with_elevs:
                ax = axin
                plot_bens_radar(stacked_data, x=total_dist / 1000.0, y=y, elev=elev, ax=ax)
                fig2, ax2 = plot_bens_radar(stacked_data, x=total_dist / 1000.0, y=y, elev=None)
            else:
                fig, ax = plot_bens_radar(stacked_data, x=total_dist / 1000.0, y=y, elev=elev)
                ax2 = axin
                plot_bens_radar(stacked_data, x=total_dist / 1000.0, y=y, elev=None, ax=ax2)
        else:
            ax = axin
            fig2, ax2 = plot_bens_radar(stacked_data, x=total_dist / 1000.0, y=y, elev=None)
    
    bb, ab = signal.butter(2, 0.05)
    ldict = {}
    if layers is not None:
        for linefn in layers:
            with codecs.open(linefn, 'r', 'UTF-8') as fin:
                lines = fin.readlines()
            la = np.array(list(map(lambda line: list(map(lambda x: np.nan if len(x) == 0 else float(x), line.rstrip('\n\r').split(','))), lines[1:])))
            indices = lines[0].rstrip('\n').split(',')

            dist = np.empty((la.shape[0], ))
            elevs = np.empty((la.shape[0], ))
            coords = np.empty((la.shape[0], 2))
            for i in range(la.shape[0]):
                dv = np.where(np.logical_and(np.abs(total_lat.flatten() - la[i, indices.index('Lat')]) < 1.0e-5, np.abs(total_lon.flatten() - la[i, indices.index('Long')]) < 1.0e-5))[0] 
                if len(dv) == 0:
                    dist[i] = np.nan
                    elevs[i] = np.nan
                    coords[i, 0] = np.nan
                    coords[i, 1] = np.nan
                else:
                    dist[i] = total_dist.flatten()[dv[0]]
                    elevs[i] = elev[dv[0]]
                    coords[i, 0] = total_lon.flatten()[dv[0]] 
                    coords[i, 1] = total_lat.flatten()[dv[0]] 
            if plot_layer is None:
                lrange = range(1, 19)
            else:
                lrange = [plot_layer]

            for linenum in lrange:
                name = 'Layer {:d} 2-Way Time'.format(linenum)
                if name not in indices:
                    continue
                else:
                    depth = tt_to_m_variable_arr(diels, la[:, indices.index(name)])
                    depth[depth == 0] = np.nan
                    try:
                        depth[~np.isnan(depth)] = signal.filtfilt(bb, ab, depth[~np.isnan(depth)])
                    except ValueError:
                        pass

                    valid_mask = np.logical_and(~np.isnan(dist), ~np.isnan(depth))
                    if np.any(valid_mask):
                        slopeval, intercept, _, _, _ = linregress(dist[valid_mask] / 1000., depth[valid_mask])
                    else:
                        slopeval, intercept = np.nan, np.nan

                    if axin is None:
                        with open('line_{:d}.csv'.format(linenum), 'w') as fout:
                            fout.write('lon, lat, distance, depth, elevation\n')
                            for i in range(len(depth)):
                                fout.write('{:f}, {:f}, {:f}, {:f}, {:f}\n'.format(coords[i, 0], coords[i, 1], dist[i], depth[i], elevs[i])) 

                    ldict['layer {:d}'.format(linenum)] = np.hstack((coords, np.vstack((dist, elevs, depth)).transpose()))
                    
                    if plot_layer is None:
                        print(dist / 1000., elevs - depth)
                        pl = ax.plot(dist / 1000., elevs - depth, linewidth=1)
                    else:
                        pl = ax.plot(dist / 1000., elevs - depth, linewidth=1, color='C0')
                    c = pl[0].get_color()
                    if slope:
                        da = np.array([0, 100])
                        ax.plot(da, slopeval * da + intercept, color=pl[0].get_color(), linestyle='-')
                        print(da, slopeval * da + intercept)
                        print('{:e} {:d}'.format(slopeval, linenum))

                    try:
                        up50_loc = (-89.539132, 137.130607)
                        udist = (la[:, indices.index('Lat')] - up50_loc[0]) ** 2.0 + (la[:, indices.index('Long')] - up50_loc[1]) ** 2.0 / 1.0e6  # this is scaled in longitude b/c we are so close to pole that it will dominate o/w
                        up50 = np.argmin(udist[valid_mask])
                    except:
                        up50 = -1

                    if axin is None:
                        if np.sum(valid_mask) > 0:
                            ax.plot(dist[valid_mask][up50] / 1000., elevs[valid_mask][up50] - depth[valid_mask][up50], linestyle='none', marker='*', color=c, markersize=markersize)
                            ax.plot(dist[valid_mask][0] / 1000., elevs[valid_mask][0] - depth[valid_mask][0], linestyle='none', marker='*', color=c, markersize=markersize)
                            print('Layer {:d} at pole has depth {:f}'.format(linenum, depth[~np.isnan(depth)][0]))
                            print('Layer {:d} at USP has depth {:f}'.format(linenum, depth[~np.isnan(depth)][up50]))

                    if bold_layer is None or linenum != bold_layer:
                        lwn = lw
                    else:
                        lwn = lw * 3
                    pl = ax2.plot(dist / 1000., depth, linewidth=lwn)

                    if label:
                        if linenum < 15:
                            ln = linenum + 3
                        else:
                            ln = linenum - 15
                        if len(depth[valid_mask]) > 0:
                            ax2.text(-1.5, depth[valid_mask][0], '{:d}'.format(ln), color=pl[0].get_color(), fontsize=8, va='center', ha='center')
                        
                    c = pl[0].get_color()
                    if slope:
                        ax2.plot(dist / 1000., dist / 1000. * slopeval + intercept, color=pl[0].get_color(), linestyle='-')

                    try:
                        up50_loc = (-89.539132, 137.130607)
                        udist = (la[:, indices.index('Lat')] - up50_loc[0]) ** 2.0 + (la[:, indices.index('Long')] - up50_loc[1]) ** 2.0 / 1.0e6  # this is scaled in longitude b/c we are so close to pole that it will dominate o/w
                        if np.min(udist[valid_mask]) < 0.001:
                            up50 = np.argmin(udist[valid_mask])
                        else:
                            up50 = 99999999
                    except:
                        up50 = -1

                    try:
                        ax2.plot(dist[valid_mask][0] / 1000., depth[valid_mask][0], linestyle='none', marker='*', color=c, markersize=markersize)
                        ax2.plot(dist[valid_mask][up50] / 1000., depth[valid_mask][up50], linestyle='none', marker='^', color=c, markersize=markersize)
                    except IndexError:
                        # if there is no data, don't worry!
                        pass

    if axin is None:
        fig.savefig('test.png', dpi=400)
        fig2.savefig(os.path.splitext(fns[0])[0] + '.pdf', dpi=400)
    
    return (gps_data, stacked_data, kinematic_data, total_lat, total_lon, total_dist, dzts, elev_list), ldict


def apply_gain(gainpoints, stacked_data):
    gp = np.atleast_2d(np.array(list(map(float, gainpoints.split(','))))).transpose()
    space = stacked_data.shape[0] // (len(gp) - 1)
    pad_size = stacked_data.shape[0] - (len(gp) - 1) * space
    gain = np.atleast_2d(np.hstack([np.linspace(gp[i], gp[i + 1], num=space) for i in range(len(gp) - 1)] + [gp[-1] for i in range(pad_size)])).transpose()
    # Now we convert the gain to decibels
    gain = 10. ** (gain / 10.)
    stacked_data *= gain
    return stacked_data


def gained_decibels(gainpoints, stacked_data, p1=None, recenter=True):
    data = data_to_db(stacked_data, p1=p1)
    data = apply_gain(gainpoints, data)
    if recenter:
        data -= data[data.shape[0] // 2:, :].mean()
    return data


def zero_surface(data):
    avg_val = np.nanmean(data, axis=1)
    max_ind = np.argmax(avg_val)
    return data[max_ind:, :], max_ind


def plot_bens_radar(data, x=None, y=None, out_fn=None, elev=None, ax=None, h_hpass=0.01, h_lpass=0.1, v_hpass=0.1, v_lpass=0.2):
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    else:
        fig = None
    X, Y = np.meshgrid(x, y)

    # now let's do some filtering, start with both high and low passes
    if h_hpass is not None:
        b, a = signal.butter(2, h_hpass, btype='high')
        data = signal.filtfilt(b, a, data, axis=1)
    if h_lpass is not None:
        b, a = signal.butter(2, h_lpass, btype='low')
        data = signal.filtfilt(b, a, data, axis=1)
    # b_vert, a_vert = signal.iirfilter(11, [0.0001, 200], rs=60, btype='band', analog=True, ftype='cheby2')

    if v_lpass is not None:
        b_hor, a_hor = signal.butter(2, v_lpass)
        data = signal.filtfilt(b_hor, a_hor, data, axis=0)
    if v_hpass is not None:
        b_hor, a_hor = signal.butter(2, v_hpass, btype='high')
        data = signal.filtfilt(b_hor, a_hor, data, axis=0)

    lims = np.percentile(data, (10, 90))

    if elev is not None:
        Y *= -1
        Y += elev
    if x is not None and y is not None:
        # plt.contourf(x, y, np.flipud(data), 2048, cmap=plt.cm.gray_r, norm=LogNorm(vmin=data[color_subset_minind:, :].min(), vmax=data[color_subset_minind:, :].max()))
        # levels = np.linspace(data.min(), data.max(), num=2048)
        try:
            ax.pcolormesh(X, Y, np.flipud(data), cmap=plt.cm.gray_r, vmin=lims[0], vmax=lims[1])
        except:
            ax.pcolormesh(X[1:], Y, np.flipud(data), cmap=plt.cm.gray_r, vmin=lims[0], vmax=lims[1])
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(Y.min(), Y.max())
        if elev is None:
            ax.invert_yaxis()
    else:
        # plt.imshow(data, cmap=plt.cm.gray_r, norm=LogNorm(vmin=data.min(), vmax=data.max()))
        # plt.imshow(data, cmap=plt.cm.bwr, norm=LogNorm(vmin=1.0e-6, vmax=data.max()))
        # plt.imshow(data, cmap=plt.cm.gray_r, vmin=data.min(), vmax=data.max())
        plt.imshow(data, cmap=plt.cm.gray_r, vmin=data.min() / 10., vmax=data.max() / 10.)
    if out_fn is not None:
        plt.savefig(out_fn)
    return fig, ax


def data_to_db(data, p1=None):
    if p1 is None:
        shp = data.shape
        p1 = data[shp[0] // 2:, :].mean()
    return 10. * np.log(data / p1)


def read_mat():
    mat = loadmat('kinematic_elevations.mat')
    return mat['x'], mat['y'], mat['elev']


class GssiError(Exception):
    pass
