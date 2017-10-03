#! /usr/bin/env python

import os
import argparse

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from ..lib.gssilib import read_dzt


def plot_bens_radar(dat, out_fn=None):
    data = dat.samp
    plt.figure(figsize=(12, 8))
    plt.imshow(data, cmap=plt.cm.gray_r, norm=LogNorm(vmin=data.min(), vmax=data.max()))
    if out_fn is not None:
        plt.savefig(out_fn)
    else:
        plt.show()


def main():
    parser = get_args()
    args = parser.parse_args()
    for fn in args.fn:
        plot_fn(fn)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fn', nargs='+', type=str, help='Files to process')
    return parser


def plot_fn(fn):
    dat = read_dzt(fn)
    plot_bens_radar(dat, out_fn=os.path.splitext(fn)[0] + '.eps')


if __name__ == '__main__':
    main()
