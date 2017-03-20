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
import os
import pickle
import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm, LinearSegmentedColormap

# import gpslib
import gssilib

colors = [(0, 0, 1), (1, 1, 1), (1, 1, 1), (1, 1, 1), (1, 0, 0)]  # B W R
bwor = LinearSegmentedColormap.from_list('wide_bwor', colors)


def main():
    parser = get_args()
    args = parser.parse_args()
    check_args(args)
    rev_list = [False if i == '0' else True for i in args.rev.split(',')]
    dzts = [gssilib.read_dzt(fn, rev=rev) for fn, rev in zip(args.fn, rev_list)]
    gssilib.check_headers(dzts)
    gps_data = [gssilib.get_dzg_data(os.path.splitext(fn)[0] + '.DZG', args.t_srs, rev=rev) for fn, rev in zip(args.fn, rev_list)]

    # detrend horizontally because of weird effects, then normalize so the different chunks match as well
    if True:
        for dzt in dzts:
            mean_vals = np.mean(dzt.samp[dzt.samp.shape[0] // 2:, :], axis=0)
            mean_val = np.mean(mean_vals)
            mean_vals /= mean_val
            vals = np.vstack((np.ones(mean_vals.shape), np.arange(len(mean_vals)))).transpose()
            fit = np.linalg.lstsq(vals, mean_vals)[0]
            dzt.samp = dzt.samp / np.dot(vals, fit)
        all_mean_vals = np.array([np.mean(dzt.samp[dzt.samp.shape[0] // 2:, :]) for dzt in dzts])
        meanest_value = np.mean(all_mean_vals)
        all_mean_vals /= meanest_value
        for dzt, scale  in zip(dzts, all_mean_vals):
            dzt.samp = dzt.samp / scale

    gps_stack_number = (gps_data[0].scans[1] - gps_data[0].scans[0]) * args.stack
    stacked_data_list =[np.array([np.nanmean(dzts[j].samp[:, i * gps_stack_number:(i + 1) * gps_stack_number], axis=1) for i in range(dzts[j].samp.shape[1] // gps_stack_number)]).transpose() for j in range(len(dzts))]

    stacked_data = np.hstack(stacked_data_list)
    stacked_data = zero_surface(stacked_data)
    if args.gp is not None:
        stacked_data = gained_decibels(args.gp, stacked_data)
    else:
        stacked_data = data_to_db(stacked_data)

    # Now find the x coordinates for plotting
    for gpd in gps_data:
        gpd.set_proj('sps')
        gpd.get_dist()

    dist_list = [gps_data[0].dist.flatten()[1:]] + [gps_data[i].dist.flatten()[1:] for i in range(1, len(gps_data))]
    dist_list = [d[::args.stack] for d in dist_list]

    # for some reason my other approach was buggy, so I'm being lazy
    for i in range(1, len(dist_list)):
        dist_list[i] += dist_list[i-1][-1]
    for i, dzt in enumerate(stacked_data_list):
        if dzt.shape[1] != len(dist_list[i]):
            dist_list[i] = dist_list[i][:-1]
    total_dist = np.hstack(dist_list)

    if args.diel is not None:
        # Find the travel time per pixel
        incremental_travel_time = dzts[0].header.range * 1.0e-9 * 3.0e8 / np.sqrt(args.diel) / 2. / dzts[0].header.nsamp
        y = np.flipud(np.arange(stacked_data.shape[0])) * incremental_travel_time
    else:
        y = np.flipud(np.arange(stacked_data.shape[0]))
    plot_bens_radar(stacked_data, x=total_dist / 1000.0, y=y, out_fn=args.o + '.png')


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
    return data[max_ind:, :]


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gp', type=str, default=None, help='gainpoints, in a comma separated list without spaces')
    parser.add_argument('--stack', type=int, default=1, help='stack this many points (already interpolated to gps points)')
    parser.add_argument('--diel', type=float, default=2.6, help='dielectric constant')
    parser.add_argument('fn', type=str, nargs='+', help='DZT filenames to use and mash together, in order')
    parser.add_argument('--rev', type=str, default=None, help='reverse files, comma separated list of 0 not to flip and 1 to flip')
    parser.add_argument('--res', type=float, default=20., help='resample the data to this spatial resolution')
    parser.add_argument('--t_srs', type=str, default='sps', choices=['sps', 'wgs84', 'EPSG:3031', 'll'], help='target projection')
    parser.add_argument('--o', type=str, default=None, help='Prefix for outputs')
    return parser


def plot_bens_radar(data, x=None, y=None, out_fn=None):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.gca()
    lims = np.percentile(data, (1, 99))
    av_lim = np.mean(np.abs(lims))
    if x is not None and y is not None:
        # plt.contourf(x, y, np.flipud(data), 2048, cmap=plt.cm.gray_r, norm=LogNorm(vmin=data[color_subset_minind:, :].min(), vmax=data[color_subset_minind:, :].max()))
        # levels = np.linspace(data.min(), data.max(), num=2048)
        # plt.pcolormesh(x, y, np.flipud(data), cmap=plt.cm.gray_r, vmin=lims[0], vmax=lims[1])
        plt.pcolormesh(x, y, np.flipud(data), cmap=bwor, norm=SymLogNorm(1.0, vmin=-av_lim, vmax=av_lim))
        ax.set_xlabel('Distance (km)')
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(y.min(), y.max())
        ax.invert_yaxis()
    else:
        # plt.imshow(data, cmap=plt.cm.gray_r, norm=LogNorm(vmin=data.min(), vmax=data.max()))
        # plt.imshow(data, cmap=plt.cm.bwr, norm=LogNorm(vmin=1.0e-6, vmax=data.max()))
        # plt.imshow(data, cmap=plt.cm.gray_r, vmin=data.min(), vmax=data.max())
        plt.imshow(data, cmap=plt.cm.gray_r, vmin=data.min() / 10, vmax=data.max() / 10)
    if out_fn is not None:
        plt.savefig(out_fn, dpi=400)
    else:
        plt.show()


def data_to_db(data, p1=None):
    if p1 is None:
        shp = data.shape
        p1 = data[shp[0] // 2:, :].mean()
    dat = 10. * np.log(data / p1)
    return 10. * np.log(data / p1)


def check_args(args):
    if args.gp is not None:
        try:
            list(map(float, args.gp.split(',')))
        except:
            raise gssilib.GssiError('Gainpoints must be a comma separated list')


if __name__ == '__main__':
    main()
