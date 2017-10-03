#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2016 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""
Combine GPS and DZT info to form a nice output
"""

import tracelib

import os
import codecs
import pickle
import argparse
import hashlib

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.stats import linregress
from scipy.io import loadmat
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Button

import gpslib
import gssilib


def main():
    parser = get_args()
    args = parser.parse_args()
    check_args(args)
    hashval = hashlib.sha256(''.join(args.fn).encode('UTF-8')).hexdigest()
    if args.pickle is not None or os.path.exists(hashval):
        if args.pickle is not None:
            hashval = args.pickle
        (gps_data,
         stacked_data,
         kinematic_data,
         total_lat,
         total_lon,
         total_dist,
         dzts,
         elev_list) = pickle.load(open(hashval, 'rb'))

    else:
        rev_list = [False if i == '0' else True for i in args.rev.split(',')]
        gps_data = [gssilib.get_dzg_data(os.path.splitext(fn)[0] + '.DZG', args.t_srs, rev=rev) for fn, rev in zip(args.fn, rev_list)]
        # Now find the x coordinates for plotting
        for gpd in gps_data:
            gpd.set_proj('sps')
            gpd.get_dist()

        if args.elev is not None:
            kinematic_data = gpslib.kinematic_info('kinematic_elevations.mat')
            elev_list = [kinematic_data.query(gpd.x, gpd.y) for gpd in gps_data]
        else:
            kinematic_data = None
            elev_list = None
        total_dist = np.hstack([gps_data[0].dist.flatten()[0:]] + [gps_data[i].dist.flatten()[1:] + np.sum([gps_data[j].dist.flatten()[-1] for j in range(i)]) for i in range(1, len(gps_data))])
        total_lat = np.hstack([gps_data[0].lat.flatten()[0:]] + [gps_data[i].lat.flatten()[1:] for i in range(1, len(gps_data))])
        total_lon = np.hstack([gps_data[0].lon.flatten()[0:]] + [gps_data[i].lon.flatten()[1:] for i in range(1, len(gps_data))])
        dzts = [gssilib.read_dzt(fn, rev=rev) for fn, rev in zip(args.fn, rev_list)]
        # dzts = pickle.load(open('pickled_dzt', 'rb'))
        gssilib.check_headers(dzts)
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
        print('dumped')

    if args.elev is not None:
        elev = np.concatenate(elev_list).flatten()[:len(total_dist)]
    else:
        elev = np.zeros(total_dist.shape).flatten()

    stacked_data, zero_ind = zero_surface(stacked_data)

    if args.gp is not None:
        stacked_data = gained_decibels(args.gp, stacked_data)
    else:
        stacked_data = data_to_db(stacked_data)

    if args.dielf is not None:
        diels = gssilib.load_dielf(args.dielf)
        
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
    elif args.diel is not None:
        # Find the travel time per pixel
        incremental_travel_time = dzts[0].header.range * 1.0e-9 * 3.0e8 / np.sqrt(args.diel) / 2. / dzts[0].header.nsamp
        y = np.flipud(np.arange(stacked_data.shape[0])) * incremental_travel_time
    else:
        y = np.flipud(np.arange(stacked_data.shape[0]))

    ax, buttons = plot_bens_radar(stacked_data, x=total_dist / 1000.0, y=y, elev=elev)
    zf = tracelib.zoom_factory(ax)
    zf.connect()

    trace_layers(ax, x=total_dist / 1000.0, y=y, elev=elev, buttons=buttons)
    
    hann = signal.hanning(21)
    if args.layers is not None:
        with open('layer.csv', 'w') as fout:
            with open('ulayer.csv', 'w') as fout2:
                for linefn in args.layers:
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

                    for linenum in range(1, 20):
                        name = 'Layer {:d} 2-Way Time'.format(linenum)
                        if name not in indices:
                            continue
                        else:
                            depth = gssilib.tt_to_m_variable_arr(diels, la[:, indices.index(name)])
                            depth[depth == 0] = np.nan
                            try:
                                depth[~np.isnan(depth)][10:-10] = signal.convolve(depth[~np.isnan(depth)], hann, mode='valid', method='direct') / np.sum(hann)
                            except ValueError:
                                pass

                            valid_mask = np.logical_and(~np.isnan(dist), ~np.isnan(depth))
                            if np.any(valid_mask):
                                slope, intercept, _, _, _ = linregress(dist[valid_mask] / 1000., depth[valid_mask])
                            else:
                                slope, intercept = np.nan, np.nan
                            with open('line_{:d}.csv'.format(linenum), 'w') as fout:
                                fout.write('lon, lat, distance, depth, elevation\n')
                                for i in range(len(depth)):
                                    fout.write('{:f}, {:f}, {:f}, {:f}, {:f}\n'.format(coords[i, 0], coords[i, 1], dist[i], depth[i], elevs[i])) 

                            pl = ax.plot(dist / 1000., elevs - depth, linewidth=1)
                            c = pl[0].get_color()
                            if args.slope:
                                ax.plot(dist / 1000., dist / 1000. * slope + intercept, color=pl[0].get_color(), linestyle='-')

                            try:
                                up50_loc = (-89.539132, 137.130607)
                                udist = (la[:, indices.index('Lat')] - up50_loc[0]) ** 2.0 + (la[:, indices.index('Long')] - up50_loc[1]) ** 2.0 / 1.0e6  # this is scaled in longitude b/c we are so close to pole that it will dominate o/w
                                up50 = np.argmin(udist[valid_mask])
                            except:
                                up50 = -1

                            try:
                                ax.plot(dist[valid_mask][up50] / 1000., elevs[valid_mask][up50] - depth[valid_mask][up50], linestyle='none', marker='*', color=c)
                                ax.plot(dist[valid_mask][0] / 1000., elevs[valid_mask][0] - depth[valid_mask][0], linestyle='none', marker='*', color=c)
                                # just pass on if we have nothing not nan
                                print('Layer {:d} at pole has depth {:f}'.format(linenum, depth[~np.isnan(depth)][0]))
                                print('Layer {:d} at USP has depth {:f}'.format(linenum, depth[~np.isnan(depth)][up50]))
                                fout.write('{:d},{:f}\n'.format(linenum, depth[~np.isnan(depth)][0]))
                                fout2.write('{:d},{:f}\n'.format(linenum, depth[~np.isnan(depth)][up50]))
                            except:
                                pass
    plt.savefig('test.png', dpi=400)


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


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gp', type=str, default=None, help='gainpoints, in a comma separated list without spaces')
    parser.add_argument('--diel', type=float, default=2.6, help='dielectric constant')
    parser.add_argument('--dielf', type=str, default=None, help='CSV file of dielectric constants, with the first column as max depth and the second as dielectric constant')
    parser.add_argument('fn', type=str, nargs='+', help='DZT filenames to use and mash together, in order')
    parser.add_argument('--rev', type=str, default=None, help='reverse files, comma separated list of 0 not to flip and 1 to flip')
    parser.add_argument('--layers', type=str, nargs='+', default=None, help='files containing layers')
    parser.add_argument('--res', type=float, default=20., help='resample the data to this spatial resolution')
    parser.add_argument('-slope', action='store_true', help='plot regression')
    parser.add_argument('--t_srs', type=str, default='sps', choices=['sps', 'wgs84', 'EPSG:3031', 'll'], help='target projection')
    parser.add_argument('--elev', type=str, default=None, help='Matlab file with kinematic data')
    parser.add_argument('--pickle', type=str, default=None, help='load this pickled data')
    return parser


def plot_bens_radar(data, x=None, y=None, out_fn=None, elev=None):
    plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(10, 2, left=0.05, bottom=0.05, right=0.99, top=0.99, wspace=0.01, hspace=0.01, width_ratios=(10, 1))
    ax = plt.subplot(gs[:, 0])
    axnew = plt.subplot(gs[0, 1])
    axnext = plt.subplot(gs[1, 1])
    axprev = plt.subplot(gs[2, 1])
    axnum = plt.subplot(gs[3, 1])
    axundo = plt.subplot(gs[4, 1])
    axsave = plt.subplot(gs[5, 1])

    # axnext = plt.axes([0.91, 0.955, 0.08, 0.04])

    bnew = Button(axnew, 'New')
    bnext = Button(axnext, 'Next')
    bundo = Button(axundo, 'Undo')
    bsave = Button(axsave, 'Save')
    bprev = Button(axprev, 'Previous')

    lims = np.percentile(data, (10, 90))
    X, Y = np.meshgrid(x, y)
    Y *= -1
    Y += elev
    if x is not None and y is not None:
        # plt.contourf(x, y, np.flipud(data), 2048, cmap=plt.cm.gray_r, norm=LogNorm(vmin=data[color_subset_minind:, :].min(), vmax=data[color_subset_minind:, :].max()))
        # levels = np.linspace(data.min(), data.max(), num=2048)
        try:
            ax.pcolormesh(X, Y, np.flipud(data), cmap=plt.cm.gray_r, vmin=lims[0], vmax=lims[1])
        except:
            ax.pcolormesh(X[1:], Y, np.flipud(data), cmap=plt.cm.gray_r, vmin=lims[0], vmax=lims[1])
        ax.set_xlabel('Distance (km)')
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(Y.min(), Y.max())
        # ax.invert_yaxis()
    else:
        plt.imshow(data, cmap=plt.cm.gray_r, norm=LogNorm(vmin=data.min(), vmax=data.max()), interpolation='nearest')
        # plt.imshow(data, cmap=plt.cm.bwr, norm=LogNorm(vmin=1.0e-6, vmax=data.max()))
        # plt.imshow(data, cmap=plt.cm.gray_r, vmin=data.min(), vmax=data.max())
        # plt.imshow(data, cmap=plt.cm.gray_r, vmin=data.min() / 10, vmax=data.max() / 10)
    if out_fn is not None:
        plt.savefig(out_fn)

    buttons = {'Next': bnext, 'Undo': bundo, 'Save': bsave, 'Previous': bprev, 'Number': axnum, 'New': bnew}
    return ax, buttons


def data_to_db(data, p1=None):
    if p1 is None:
        shp = data.shape
        p1 = data[shp[0] // 2:, :].mean()
    return 10. * np.log(data / p1)


def check_args(args):
    if args.gp is not None:
        try:
            list(map(float, args.gp.split(',')))
        except:
            raise gssilib.GssiError('Gainpoints must be a comma separated list')


def trace_layers(ax, x=None, y=None, elev=None, buttons=None):
    lines = tracelib.LineList(ax)
    if buttons is not None:
        lines.add_buttons(buttons)
    plt.show()


def read_mat():
    mat = loadmat('kinematic_elevations.mat')
    return mat['x'], mat['y'], mat['elev']


if __name__ == '__main__':
    main()
