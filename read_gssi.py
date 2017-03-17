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
from matplotlib.colors import LogNorm

# import gpslib
import gssilib


def fake_main():
    parser = get_args()
    args = parser.parse_args()
    check_args(args)
    rev_list = [False if i == '0' else True for i in args.rev.split(',')]
    dzts = [gssilib.read_dzt(fn, rev=rev) for fn, rev in zip(args.fn, rev_list)]
    gssilib.check_headers(dzts)
    pickle.dump(dzts, open('pickled_dzt', 'wb'))
    gps_data = [gssilib.get_dzg_data(os.path.splitext(fn)[0] + '.DZG', args.t_srs, rev=rev) for fn, rev in zip(args.fn, rev_list)]
    pickle.dump(gps_data, open('pickled_gps', 'wb'))


def main():
    parser = get_args()
    args = parser.parse_args()
    check_args(args)
    dzts = pickle.load(open('pickled_dzt', 'rb'))
    gssilib.check_headers(dzts)
    gps_data = pickle.load(open('pickled_gps', 'rb'))
    gps_stack_number = (gps_data[0].scans[1] - gps_data[0].scans[0]) * args.stack
    stacked_data = np.hstack([np.array([np.nanmean(dzts[j].samp[:, i * gps_stack_number:(i + 1) * gps_stack_number], axis=1) for i in range(dzts[j].samp.shape[1] // gps_stack_number)]).transpose() for j in range(len(dzts))])

    stacked_data = zero_surface(stacked_data)

    if args.gp is not None:
        stacked_data = gained_decibels(args.gp, stacked_data)
    else:
        stacked_data = data_to_db(stacked_data)

    # Now find the x coordinates for plotting
    for gpd in gps_data:
        gpd.set_proj('sps')
        gpd.get_dist()
    total_dist = np.hstack([gps_data[0].dist.flatten()[1:]] + [gps_data[i].dist.flatten()[1:] + gps_data[i - 1].dist.flatten()[-1] for i in range(1, len(gps_data))])
    total_dist = total_dist[::args.stack]

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
    lims = np.percentile(data, (5, 95))
    if x is not None and y is not None:
        # plt.contourf(x, y, np.flipud(data), 2048, cmap=plt.cm.gray_r, norm=LogNorm(vmin=data[color_subset_minind:, :].min(), vmax=data[color_subset_minind:, :].max()))
        # levels = np.linspace(data.min(), data.max(), num=2048)
        plt.pcolormesh(x, y, np.flipud(data), cmap=plt.cm.gray_r, vmin=lims[0], vmax=lims[1])
        ax.set_xlabel('Distance (km)')
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(y.min(), y.max())
        ax.invert_yaxis()
    else:
        print(data.min(), data.max())
        print(data)
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
        print(p1)
    dat = 10. * np.log(data / p1)
    print(dat[3 * shp[0] // 4:, :].mean())
    return 10. * np.log(data / p1)


def check_args(args):
    if args.gp is not None:
        try:
            list(map(float, args.gp.split(',')))
        except:
            raise gssilib.GssiError('Gainpoints must be a comma separated list')


if __name__ == '__main__':
    main()
