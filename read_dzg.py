#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2016 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""
Use this to read gps information from a dzg file
"""
import os
import numpy as np
import argparse
from modeltools.lib.glib import dict2shp_pts

from gssilib import get_dzg_data


def main():
    parser = get_args()
    args = parser.parse_args()
    basenames = [os.path.splitext(os.path.split(fn)[1])[0] for fn in args.fn]
    basename_dict = {bn: None for bn in basenames}
    for fn, bn in zip(args.fn, basenames):
        basename_dict[bn] = get_dzg_data(fn, t_srs=args.t_srs)
        dict2shp_pts(os.path.splitext(fn)[0] + '.shp', {'x': basename_dict[bn].x, 'y': basename_dict[bn].y, 'z': basename_dict[bn].z, 'scan': basename_dict[bn].scans, 'filename': [bn for i in range(len(basename_dict[bn].x))]})
    if len(args.fn) > 1:
        if args.o is None:
            args.o = os.path.split(args.fn[0])[0] + '/' + 'combined_track.shp'
        dict2shp_pts(args.o, {'x': np.hstack([basename_dict[bn].x for bn in basenames]), 'y': np.hstack([basename_dict[bn].y for bn in basenames]), 'z': np.hstack([basename_dict[bn].z for bn in basenames]), 'scans': np.hstack([basename_dict[bn].scans for bn in basenames]), 'filename': list(sum([[bn for i in range(len(basename_dict[bn].x))] for bn in basenames], []))})
    return basename_dict


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fn', type=str, nargs='+', help='files with extension dzg to process')
    parser.add_argument('-t_srs', type=str, default='sps', choices=['sps', 'wgs84', 'll', 'EPSG:3031'], help='Convert to this reference')
    parser.add_argument('-o', type=str, default=None, help='Name for combined output')
    return parser

if __name__ == '__main__':
    main()
