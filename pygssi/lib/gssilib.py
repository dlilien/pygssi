#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2016 dlilien <dal22@uw.edu>
#
# Distributed under terms of the MIT license.

"""
Class definitions and helper functions for gssi info
"""

import numpy as np
from .gpslib import kinematic_info
from .gssi_filelib import get_dzg_data, read_dzt, check_headers
from .conversionlib import gained_decibels, data_to_db, tt_to_m_variable_arr
import os
import codecs
import pickle
import hashlib

import matplotlib.pyplot as plt
from cycler import cycler
from scipy import signal

# No gray
plt.rc('axes', prop_cycle=(cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#17becf']) + cycler('linestyle', ['solid'] * 9)))


def process_radar(fns,
                  rev_list=None,
                  t_srs='sps',
                  elev_fn=None,
                  pickle_fn=None,
                  axin=None,
                  gp=None,
                  dielf=None,
                  diel=None,
                  layers=None,
                  label=False,
                  lw=1,
                  plotdata=True,
                  with_elevs=True,
                  xoff=0.0,
                  plot_layer=None,
                  scale_x=1000.0):
    """Read in raw (or pickled) radar data and do a bunch of processing on it
    
    Parameters
    ----------
    fns: str or iterable of strs
        The file(s) to process
    rev_list: list of bools, optional
        For each file, do you want to reverse it. Default is no reversal (can also handle true or false to do the same to all files)

    """
    if type(fns) is str:
        fns = [fns]
    hashval = 'radar_' + hashlib.sha256(''.join(fns).encode('UTF-8')).hexdigest()

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
        if rev_list is None:
            rev_list = [False for fn in fns]
        elif rev_list in [True, False]:
            rev_list = [rev_list for fn in fns]

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

    def zero_surface(data):
        avg_val = np.nanmean(data, axis=1)
        max_ind = np.argmax(avg_val)
        return data[max_ind:, :], max_ind

    stacked_data, zero_ind = zero_surface(stacked_data)

    if gp is not None:
        stacked_data = gained_decibels(gp, stacked_data)
    else:
        stacked_data = data_to_db(stacked_data)

    if dielf is not None:
        def load_dielf(fn):
            with open(fn) as fin:
                lines = fin.readlines()

            # pad this with some values so we don't go out of bounds
            lines = list(map(lambda x: list(map(float, x.replace('\ufeff', '').rstrip('\n\r').split(','))), lines))
            diels = np.array([[0.0, lines[0][1]]] + lines + [[np.inf, lines[-1][1]]])
            return diels
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
        fig, ax = plot_radar(stacked_data, x=total_dist / scale_x, y=y, elev=elev)
        fig2, ax2 = plot_radar(stacked_data, x=total_dist / scale_x, y=y, elev=None)
    else:
        if plotdata:
            if with_elevs:
                ax = axin
                plot_radar(stacked_data, x=total_dist / scale_x, y=y, elev=elev, ax=ax)
                fig2, ax2 = plot_radar(stacked_data, x=total_dist / scale_x, y=y, elev=None)
            else:
                fig, ax = plot_radar(stacked_data, x=total_dist / scale_x, y=y, elev=elev)
                ax2 = axin
                plot_radar(stacked_data, x=total_dist / scale_x, y=y, elev=None, ax=ax2)
        else:
            ax = axin
            fig2, ax2 = plot_radar(stacked_data, x=total_dist / scale_x, y=y, elev=None)
    
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

                    if axin is None:
                        with open('line_{:d}.csv'.format(linenum), 'w') as fout:
                            fout.write('lon, lat, distance, depth, elevation\n')
                            for i in range(len(depth)):
                                fout.write('{:f}, {:f}, {:f}, {:f}, {:f}\n'.format(coords[i, 0], coords[i, 1], dist[i], depth[i], elevs[i])) 

                    ldict['layer {:d}'.format(linenum)] = np.hstack((coords, np.vstack((dist, elevs, depth)).transpose()))
                    
                    if plot_layer is None:
                        pl = ax.plot(dist / scale_x, elevs - depth, linewidth=1)
                    else:
                        pl = ax.plot(dist / scale_x, elevs - depth, linewidth=1)

                    pl = ax2.plot(dist / scale_x, depth, linewidth=lw)

                    ax2.text(-1.5, depth[valid_mask][0], '{:d}'.format(linenum), color=pl[0].get_color(), fontsize=8, va='center', ha='center')

    if axin is None:
        fig.savefig('test.png', dpi=400)
        fig2.savefig('flat_surf.png', dpi=400)
    
    return (gps_data, stacked_data, kinematic_data, total_lat, total_lon, total_dist, dzts, elev_list), ldict


def plot_radar(data, x=None, y=None, out_fn=None, elev=None, ax=None, h_hpass=0.01, h_lpass=0.1, v_hpass=0.1, v_lpass=0.2):
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
        try:
            ax.pcolormesh(X, Y, np.flipud(data), cmap=plt.cm.gray_r, vmin=lims[0], vmax=lims[1])
        except ValueError:
            ax.pcolormesh(X[1:], Y, np.flipud(data), cmap=plt.cm.gray_r, vmin=lims[0], vmax=lims[1])
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(Y.min(), Y.max())
        if elev is None:
            ax.invert_yaxis()
    else:
        plt.imshow(data, cmap=plt.cm.gray_r, vmin=lims[0], vmax=lims[1])

    if out_fn is not None:
        plt.savefig(out_fn)
    return fig, ax
