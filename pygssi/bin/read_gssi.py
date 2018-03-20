#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2016 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""
Combine GPS and DZT info to form a nice output
"""

import argparse

from pygssi.lib import gssilib


def main():
    parser = get_args()
    args = parser.parse_args()
    check_args(args)
    rev_list = [False if i == '0' else True for i in args.rev.split(',')]
    gssilib.process_and_plot(args.fn, rev_list, t_srs=args.t_srs, elev_fn=args.elev, pickle_fn=args.pickle, axin=None, gp=args.gp, layers=args.layers, diel=args.diel)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gp', type=str, default=None, help='gainpoints, in a comma separated list without spaces')
    parser.add_argument('--diel', type=str, default=None, help='CSV file of dielectric constants, with the first column as max depth and the second as dielectric constant')
    parser.add_argument('fn', type=str, nargs='+', help='DZT filenames to use and mash together, in order')
    parser.add_argument('--rev', type=str, default=None, help='reverse files, comma separated list of 0 not to flip and 1 to flip')
    parser.add_argument('--layers', type=str, nargs='+', default=None, help='files containing layers')
    parser.add_argument('--res', type=float, default=20., help='resample the data to this spatial resolution')
    parser.add_argument('-filter', action='store_true', help='High and low pass in the vertical and horizontal')
    parser.add_argument('--t_srs', type=str, default='sps', choices=['sps', 'wgs84', 'EPSG:3031', 'nps', 'll'], help='target projection')
    parser.add_argument('--elev', type=str, default=None, help='Matlab file with kinematic data')
    parser.add_argument('--pickle', type=str, default=None, help='load this pickled data')
    return parser


def check_args(args):
    if args.gp is not None:
        try:
            list(map(float, args.gp.split(',')))
        except TypeError:
            raise Exception('Gainpoints must be a comma separated list')

    try:
        args.diel = float(args.diel)
    except ValueError:
        pass


if __name__ == '__main__':
    main()
