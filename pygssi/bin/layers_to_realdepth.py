#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2017 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""

"""
from pygeotools.lib import geolib
from pygssi.lib import gssilib
import argparse
import codecs
import numpy as np

from scipy import signal
from scipy.io import savemat


def _get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('dielf', type=str)
    parser.add_argument('fn', type=str, nargs='+')
    parser.add_argument('--t_srs', type=str, default='sps')
    parser.add_argument('--of', type=str, default='depth_corrected_layers/')
    return parser


def main():
    parser = _get_args()
    args = parser.parse_args()
    diels = gssilib.load_dielf(args.dielf)
    hann = signal.hanning(21)
    print(diels)

    matlab_dict = {}

    for linefn in args.fn:
        with codecs.open(linefn, 'r', 'UTF-8') as fin:
            lines = fin.readlines()
        la = np.array(list(map(lambda line: list(map(lambda x: np.nan if len(x) == 0 else float(x), line.rstrip('\n\r').split(','))), lines[1:])))
        indices = lines[0].rstrip('\n').split(',')

        coords = np.empty((la.shape[0], 2))
        sps_vals = geolib.ll2sps(la[:, 2], la[:, 1])
        coords = np.vstack((sps_vals[0], sps_vals[1])).transpose()

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

                with open(args.of + '/line_{:d}.csv'.format(linenum), 'w') as fout:
                    fout.write('psx, psy, depth\n')
                    for i in range(len(depth)):
                        fout.write('{:f}, {:f}, {:f}\n'.format(coords[i, 0], coords[i, 1], depth[i]))
            if linenum == 1:
                matlab_dict['psx_layers'] = coords[:, 0]
                matlab_dict['psy_layers'] = coords[:, 1]
            # else:
            #     d = np.zeros(matlab_dict['psx_layers'].shape)
            #     d[:] = np.nan
            #     for i in range(len(coords[:, 0])):
            #         d[np.logical_and(np.abs(matlab_dict['psx_layers'] - coords[i, 0]) < 0.01, np.abs(matlab_dict['psy_layers'] - coords[i, 1]) < 0.01)] = depth[i]
            #     depth = d
            if len(depth) < len(matlab_dict['psx_layers']):
                d = np.zeros(matlab_dict['psx_layers'].shape)
                d[:] = np.nan
                d[:len(depth)] = depth
            else:
                d = depth
                    
            matlab_dict['layer_{:d}'.format(linenum)] = d

        savemat(args.of + '/layers.mat', matlab_dict)


if __name__ == '__main__':
    main()
