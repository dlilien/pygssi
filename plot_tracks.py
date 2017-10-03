#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2017 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""

"""
import codecs
import numpy as np
import matplotlib.pyplot as plt
import gssilib
from modeltools.lib import basemaplib
from pygeotools.lib import geolib

height = 6

pole_kwargs = {'projection': 'stere', 'lat_0': -90, 'lon_0': 0, 'lat_ts': -71, 'llcrnrlat': -89.0, 'urcrnrlat': -89.0, 'llcrnrlon': -175., 'urcrnrlon': 75}

spice_loc = (-1.0755393e+03, -1.5422120e+02)

camp_loc = geolib.ll2sps(137.130607, -89.539132)


class PoleBasemap(basemaplib.AntarcticaBasemap):
    background_fn = '/Users/dlilien/work/antarctica_general/moa1000_2004_hp1_v1.1.tif'
    background_res = None
    _background_ds = None

    parallels = np.atleast_1d(np.arange(-90, -88.5, 0.5))
    meridians = np.atleast_1d(np.arange(0, 375, 1))

    def __init__(self, *args, **kwargs):
        super(PoleBasemap, self).__init__(*args, **dict(pole_kwargs, **kwargs))

    def label_ll(self, plabels=[True, False, False, False], mlabels=[True, True, False, True], xoff=-0.03, yoff=-0.025, latmax=80., color='grey', linecolor='grey', mer_side_offset=-0.02, backgroundcolor='none', *args, **kwargs):
        super().label_ll(plabels=plabels, mlabels=mlabels, xoff=xoff, yoff=yoff, latmax=latmax, color=color, linecolor=linecolor, mer_side_offset=mer_side_offset, backgroundcolor=backgroundcolor, *args, **kwargs)


stake_loc_fn = '../GPS_survey/SouthPoleGISData/SP_stake_locations_xy.txt'
with open(stake_loc_fn) as fin:
    stake_locs = np.array(list(map(lambda x: list(map(float, x.rstrip('\n').lstrip(' ').split('  '))), fin.readlines())))


def plot_locs(bm):
    bm.projected_plot(stake_locs[:, 0], stake_locs[:, 1], linestyle='none', marker='.', color='k')
    bm.projected_plot(*spice_loc, marker='*')
    bm.projected_plot(camp_loc[0], camp_loc[1], marker='^')


plt.figure(figsize=(6, 6))
bm = PoleBasemap()
bm.plot_shpfile('shps/SOUTHPOLE__002.shp')
bm.plot_shpfile('shps/SOUTHPOLE__003.shp')
bm.plot_shpfile('shps/SOUTHPOLE__004.shp')
bm.drawparallels(bm.parallels, latmax=89.9)
plt.savefig('tracks.png', dpi=300)

asp = (bm.urcrnrx - bm.llcrnrx) / (bm.urcrnry - bm.llcrnry)

diels = gssilib.load_dielf('dielectric_constants_20170710.csv')

linefn = 'shallowfull.csv'
with codecs.open(linefn, 'r', 'UTF-8') as fin:
    lines = fin.readlines()
la = np.array(list(map(lambda line: list(map(lambda x: np.nan if len(x) == 0 else float(x), line.rstrip('\n\r').split(','))), lines[1:])))
indices = lines[0].rstrip('\n').split(',')
zs = {}
for i in range(1, 8):
    name = 'Layer {:d} 2-Way Time'.format(i)
    z = gssilib.tt_to_m_variable_arr(diels, la[:, indices.index(name)])
    z[z == 0] = np.nan
    zs[name] = z
x, y, z = geolib.ll2sps(la[:, indices.index('Long')], la[:, indices.index('Lat')])

linefn = 'deepfull.csv'
with codecs.open(linefn, 'r', 'UTF-8') as fin:
    lines = fin.readlines()
la2 = np.array(list(map(lambda line: list(map(lambda x: np.nan if len(x) == 0 else float(x), line.rstrip('\n\r').split(','))), lines[1:])))
indices2 = lines[0].rstrip('\n').split(',')
for i in range(8, 15):
    name = 'Layer {:d} 2-Way Time'.format(i)
    z = gssilib.tt_to_m_variable_arr(diels, la2[:, indices2.index(name)])
    z[z == 0] = np.nan
    zs[name] = z
x2, y2, z2 = geolib.ll2sps(la2[:, indices.index('Long')], la2[:, indices.index('Lat')])

# Now let's do some binned stuff
zpercs = {name: z / np.nanmean(z) for name, z in zs.items()}
zdemeaned = {name: z - np.nanmean(z) for name, z in zs.items()}
zpole = {name: z[~np.isnan(z)][0] - z for name, z in zs.items()}
zperpole = {name: z / z[~np.isnan(z)][0] for name, z in zs.items()}

if False:
    for i in range(1, 15):
        fig = plt.figure(figsize=(height * asp, height))
        ax = fig.gca()
        bm = PoleBasemap(ax=ax)
        name = 'Layer {:d} 2-Way Time'.format(i)
        cm = bm.projected_scatter(x, y, 2, c=zs[name])
        bm.fancy_colorbar(cm)
        fig.tight_layout(pad=0)
        fig.savefig('layerdepth_{:d}.png'.format(i), dpi=400)

if True:
    for name, zp in zpole.items():
        fig = plt.figure(figsize=(height * asp, height))
        i = int(name.split(' ')[1])
        ax = fig.gca()
        bm = PoleBasemap(ax=ax)
        plot_locs(bm)
        cm = bm.projected_scatter(x, y, 2, c=zp, vmin=-5, vmax=5, cmap=plt.get_cmap('BrBG'))
        bbr = bm.fancy_colorbar(cm)
        bbr.set_label('Mean layer change from pole (m)')
        fig.tight_layout(pad=0)
        fig.savefig('change_pole_{:d}.png'.format(i), dpi=400)

# overall
fig = plt.figure(figsize=(height * asp, height))
ax = fig.gca()
bm = PoleBasemap(ax=ax)
bm.moa()
plot_locs(bm)
cm = bm.projected_scatter(x, y, 2, c=np.nanmean(np.array([z for z in zpercs.values()]), axis=0) * 100 - 100, cmap=plt.get_cmap('BrBG'), vmin=-5, vmax=5)
cbr = bm.fancy_colorbar(cm, extend='both')
cbr.set_label(r'\% difference from mean')
fig.tight_layout(pad=0)
fig.savefig('perchange_mean.png', dpi=400)

fig = plt.figure(figsize=(height * asp, height))
ax = fig.gca()
bm = PoleBasemap(ax=ax)
plot_locs(bm)
cm = bm.projected_scatter(x, y, 2, c=np.nanmean(np.array([z for z in zdemeaned.values()]), axis=0))
bm.fancy_colorbar(cm)
fig.tight_layout(pad=0)
fig.savefig('change_demean.png', dpi=400)

fig = plt.figure(figsize=(height * asp, height))
ax = fig.gca()
bm = PoleBasemap(ax=ax)
plot_locs(bm)
cm = bm.projected_scatter(x, y, 2, c=np.nanmean(np.array([z for z in zpole.values()]), axis=0), vmin=-10, vmax=10, cmap=plt.get_cmap('BrBG'))
bbr = bm.fancy_colorbar(cm)
bbr.set_label('Mean layer change from pole (m)')
fig.tight_layout(pad=0)
fig.savefig('change_pole.png', dpi=400)

fig = plt.figure(figsize=(height * asp, height))
ax = fig.gca()
bm = PoleBasemap(ax=ax)
plot_locs(bm)
cm = bm.projected_scatter(x, y, 2, c=np.nanmean(np.array([z for name, z in zperpole.items() if name in indices2]), axis=0) * 100 - 100., cmap=plt.get_cmap('BrBG'), vmin=-10, vmax=10)
cbr = bm.fancy_colorbar(cm, extend='both')
cbr.set_label(r'Mean layer depth change from SPICE (\%)')
fig.tight_layout(pad=0)
fig.savefig('changeper_pole.png', dpi=400)
