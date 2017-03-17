#! /usr/bin/env python3
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from gssilib import read_dzt


def plot_bens_radar(dat, out_fn=None):
    data = dat.samp
    plt.figure(figsize=(12, 8))
    plt.imshow(data, cmap=plt.cm.gray_r, norm=LogNorm(vmin=data.min(), vmax=data.max()))
    if out_fn is not None:
        plt.savefig(out_fn)
    else:
        plt.show()


def main(fn):
    dat = read_dzt(fn)
    plot_bens_radar(dat, out_fn=os.path.splitext(fn)[0] + '.eps')


if __name__ == '__main__':
    for fn in sys.argv[1:]:
        main(fn)
