#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2016 dlilien <dal22@uw.edu>
#
# Distributed under terms of the MIT license.

"""
Class definitions and helper functions for gssi info

The workhorses here are :func:`~pygssi.lib.gssilib.process_and_plot` and :func:`~pygssi.lib.gssilib.process_radar`, which will handle a lot of the steps of taking a raw DZG file and plotting it and/or returning it in a useful form
"""

import numpy as np
from .conversionlib import gained_decibels, data_to_db, tt_to_m_variable_arr
from .gssi_filelib import read
import codecs

import matplotlib.pyplot as plt
from cycler import cycler
from scipy import signal

# No gray
plt.rc('axes', prop_cycle=(cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#17becf']) + cycler('linestyle', ['solid'] * 9)))


def process_and_plot(fns,
                     rev_list=None,
                     t_srs='sps',
                     elev_fn=None,
                     pickle_fn=None,
                     gp=None,
                     diel=None,
                     xoff=0.0,
                     filter_dat=True,
                     axin=None,
                     with_elevs=True,
                     layers=None,
                     lw=1,
                     scale_x=1000.0):
    """Read in raw (or pickled) radar data and do a bunch of processing on it, then plot it up

    This will always do two plots, one with a variable surface and one without, though if there is no elevation information the one with a variable surface is not saved.
    
    Parameters
    ----------
    fns: str or iterable of strs
        The file(s) to process
    rev_list: list of bools, optional
        For each file, do you want to reverse it. Default is no reversal (can also handle true or false to do the same to all files)
    t_srs: str, optional
        Try to project the output to this coordinate system
    elev_fn: str, optional
        A .mat file containing fields x, y, z for the surface of the radar
    gp: str, optional
        A string of comma-separated floats for the gainpoints. These get evenly distributed across depths. They are applied exponentially
    diel: float, numpy.ndarray, str, or None
        diel can be a 2-column csv with depths and dielectric constants, an array with depths and dielectric constants, a number giving a depth-constant dielectric constant, or None, which leaves things in two-way travel time. If layer_time_increment is None, we always stay in TWTT
    xoff: float, optional
        offset all the x values by this much (useful for plotting)
    filter_dat: bool, optional
        hit the data with the defaults of :func:`~pygssi.lib.gssilib.filter_data`. If you want custom values, you should take the data returned and hit it with the custom filter.
    axin: matplotlib.pyplot.Axis, optional
        An axis to plot on. If with_elevs, then this is the variable surface axis. If not, it is the flat axis. 
    with_elevs: bool, optional
        If True (default), the axis passed in is for the variable surface.
    layers: list, optional
        List of files containing layers to plot. This is designed to be used in conjunction with the radan tracing software, so the format is a 1-line header with "Dist.(m),Lat,Long,Elev(m),Time,Layer 1 2-Way Time,Layer 2 2-Way Time...".
    lw: float, optional
        Line width of layers plotted.
    scale_x: float, optional
        Divide distance by this amount for x axis. Default is 1000.0 to make polar stereographic plots in km.


    Returns
    -------
    data: nxm numpy.ndarray
        The radar data (in decibels, possibly gained)
    lldist: mx3 numpy.ndarray
        The lon, lat, and distance for each point. If no GPS data is given, each row is 0, 1, 2...
    elev: mx1 numpy.ndarray
        Elevations for each trace if kinematic GPS data was input, else None
    z: nx1 numpy.ndarray
        The depth of each point. Just 0, 1... if no dielectric constant is given
    diels: lx2 numpy.ndarray
        Depth, dielectric constant pairs.
    ldict: dict
        A dictionary of layer names and depths. None if no layers are input
    """
    data, lldist, elev, y, diels = process_radar(fns, rev_list=rev_list, t_srs=t_srs, pickle_fn=pickle_fn, elev_fn=elev_fn, gp=gp, diel=diel, xoff=xoff, filter_dat=filter_dat)
    if axin is None:
        fig, ax = plot_radar(data, x=lldist[:, 2] / scale_x, y=y, elev=elev)
        fig2, ax2 = plot_radar(data, x=lldist[:, 2] / scale_x, y=y, elev=None)
    else:
        if with_elevs:
            ax = axin
            plot_radar(data, x=lldist[:, 2] / scale_x, y=y, elev=elev, ax=ax)
            fig2, ax2 = plot_radar(data, x=lldist[:, 2] / scale_x, y=y, elev=None)
        else:
            fig, ax = plot_radar(data, x=lldist[:, 2] / scale_x, y=y, elev=elev)
            ax2 = axin
            plot_radar(data, x=lldist[:, 2] / scale_x, y=y, elev=None, ax=ax2)

    if layers is not None:
        ldict = plot_layers(layers, ax2, lldist, diels, ax_variablesurf=ax, linewidth=lw, scale_x=scale_x)
    else:
        ldict = None

    if axin is None:
        if elev is not None:
            fig.savefig('radar_variable_surf.png', dpi=400)
        fig2.savefig('radar_flat_surf.png', dpi=400)

    return data, lldist, elev, y, diels, ldict


def process_radar(fns,
                  rev_list=None,
                  t_srs='sps',
                  elev_fn=None,
                  pickle_fn=None,
                  gp=None,
                  diel=None,
                  xoff=0.0,
                  filter_dat=True,
                  cache=False):
    """Read in raw (or pickled) radar data and do a bunch of processing on it
    
    Parameters
    ----------
    fns: str or iterable of strs
        The file(s) to process
    rev_list: list of bools, optional
        For each file, do you want to reverse it. Default is no reversal (can also handle true or false to do the same to all files)
    t_srs: str, optional
        Try to project the output to this coordinate system
    elev_fn: str, optional
        A .mat file containing fields x, y, z for the surface of the radar
    gp: str, optional
        A string of comma-separated floats for the gainpoints. These get evenly distributed across depths. They are applied exponentially
    diel: float, numpy.ndarray, str, or None
        diel can be a 2-column csv with depths and dielectric constants, an array with depths and dielectric constants, a number giving a depth-constant dielectric constant, or None, which leaves things in two-way travel time. If layer_time_increment is None, we always stay in TWTT
    xoff: float, optional
        offset all the x values by this much (useful for plotting)
    filter_dat: bool, optional
        hit the data with the defaults of :func:`~pygssi.lib.gssilib.filter_data`. If you want custom values, you should take the data returned and hit it with the custom filter.
    cache: bool, optional
        pickle results for fast future use (default False)

    Returns
    -------
    data: nxm numpy.ndarray
        The radar data (in decibels, possibly gained)
    lldist: mx3 numpy.ndarray
        The lon, lat, and distance for each point. If no GPS data is given, each row is 0, 1, 2...
    elev: mx1 numpy.ndarray
        Elevations for each trace if kinematic GPS data was input, else None
    z: nx1 numpy.ndarray
        The depth of each point. Just 0, 1... if no dielectric constant is given
    diels: lx2 numpy.ndarray
        Depth, dielectric constant pairs.
    """
    if type(fns) is str:
        fns = [fns]
    if rev_list is None:
        rev_list = [False for fn in fns]
    elif rev_list in [True, False]:
        rev_list = [rev_list for fn in fns]

    _, data, lldist, dzts, elev_list = read(fns, rev_list, elev_fn=elev_fn, pickle_fn=pickle_fn, t_srs=t_srs, cache=cache)

    lldist[2, :] = lldist[2, :] + xoff

    if elev_list is not None:
        elev = np.concatenate(elev_list).flatten()[:lldist.shape[0]]
    else:
        elev = None

    data, zero_ind = zero_surface(data)

    if gp is not None:
        data = gained_decibels(gp, data)
    else:
        data = data_to_db(data)

    layer_time_increment = dzts[0].header.range / dzts[0].header.nsamp * 1.0e-9
    y, diels = apply_diel(diel, data, layer_time_increment=layer_time_increment)

    if filter_dat:
        data = filter_data(data)
    
    return data, lldist, elev, y, diels


def plot_layers(layer_fns, ax_flat, lldist, diels, ax_variablesurf=None, linewidth=1, scale_x=1000., z=None):
    """Plot layers on axes given

    Parameters
    ----------
    layer_fns: list
        List of files containing layers to plot. This is designed to be used in conjunction with the radan tracing software, so the format is a 1-line header with "Dist.(m),Lat,Long,Elev(m),Time,Layer 1 2-Way Time,Layer 2 2-Way Time...".
    ax_flat: matplotlib.pyplot.Axis
        Axis in depth (or TWTT) space on which to plot
    lldist: mx3 numpy.ndarray
        The lon, lat, and distance for each point. If no GPS data is given, each row is 0, 1, 2...
    diels: lx2 numpy.ndarray
        Depth, dielectric constant pairs.
    ax_variablesurf: matplotlib.pyplot.Axis, optional
        Axis in elevation space on which to plot
    linewidth: float, optional
        Line width of layers plotted.
    scale_x: float, optional
        Divide distance by this amount for x axis. Default is 1000.0 to make polar stereographic plots in km.
    z: mx1 numpy.ndarray, optional
        Elevations for each trace. Ignored if no ax_variablesurf.

    Returns
    -------
    ldict: dict
        A dictionary of layer names and depths. None if no layers are input
    """
    bb, ab = signal.butter(2, 0.05)
    ldict = {}
    for linefn in layer_fns:
        with codecs.open(linefn, 'r', 'UTF-8') as fin:
            lines = fin.readlines()
        la = np.array(list(map(lambda line: list(map(lambda x: np.nan if len(x) == 0 else float(x), line.rstrip('\n\r').split(','))), lines[1:])))
        indices = lines[0].rstrip('\n').split(',')

        dist = np.empty((la.shape[0], ))
        elevs = np.empty((la.shape[0], ))
        coords = np.empty((la.shape[0], 2))
        for i in range(la.shape[0]):
            dv = np.where(np.logical_and(np.abs(lldist[:, 1].flatten() - la[i, indices.index('Lat')]) < 1.0e-5, np.abs(lldist[:, 0].flatten() - la[i, indices.index('Long')]) < 1.0e-5))[0] 
            if len(dv) == 0:
                dist[i] = np.nan
                elevs[i] = np.nan
                coords[i, 0] = np.nan
                coords[i, 1] = np.nan
            else:
                dist[i] = dist.flatten()[dv[0]]
                elevs[i] = z[dv[0]]
                coords[i, 0] = lldist[:, 0].flatten()[dv[0]] 
                coords[i, 1] = lldist[:, 1].flatten()[dv[0]] 
        for linenum in range(50):
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

                ldict['layer {:d}'.format(linenum)] = np.hstack((coords, np.vstack((dist, elevs, depth)).transpose()))
                
                if ax_variablesurf is not None:
                    ax_variablesurf.plot(dist / scale_x, elevs - depth, linewidth=1)
                ax_flat.plot(dist / scale_x, depth, linewidth=linewidth)
    return ldict


def plot_radar(data, x=None, y=None, out_fn=None, elev=None, ax=None, xlabel='Distance (km)', ylabel='Depth (m)'):
    """Plot an array of radar data

    Parameters
    ----------
    data: numpy.ndarray
        An nxm array of values to plot
    x: numpy.ndarray, optional
        An mx1 array x values of the data to plot. If x or y is None, the x-axis will be in number of traces
    y: numpy.ndarray, optional
        An nx1 array y values of the data to plot. If x or y is None, the y-axis will be the return bin
    out_fn: str, optional
        Save to this file (default None/don't save)
    elev: numpy.ndarray
        An nx1 array y values of surface height. Depths are subtracted from this before plotting if not None. Ignored if x or y is None.
    ax: matplotlib.pyplot.Axis
        Plot on these axes rather than creating a new figure (default is None/create figure)
    xlabel: str, optional
        Labels for the x-axis (overridden if x or y are None).
    ylabel: str, optional
        Labels for the y-axis (overridden if x or y are None).


    Returns
    -------
    fig: matplotlib.pyplot.Figure
        The figure canvas (unless an axis was passed in, in which case figure is None)
    ax: matplotlib.pyplot.Axis
        The axis instance created (or that passed in returned back if one is specified
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    else:
        fig = None
    X, Y = np.meshgrid(x, y)

    # Optional variable surface
    if elev is not None:
        Y *= -1
        Y += elev

    lims = np.percentile(data, (10, 90))

    # we need to make sure that y is not just the trace number. This is a little risky if the dielectric constant is bad news, but it should be fine
    if x is not None and y is not None and not np.all(y == np.array(list(range(len(y) - 1, -1, -1)))):
        try:
            ax.pcolormesh(X, Y, np.flipud(data), cmap=plt.cm.gray_r, vmin=lims[0], vmax=lims[1])
        except ValueError:
            ax.pcolormesh(X[1:], Y, np.flipud(data), cmap=plt.cm.gray_r, vmin=lims[0], vmax=lims[1])
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(Y.min(), Y.max())
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        if elev is None:
            ax.invert_yaxis()
    else:
        plt.imshow(data, cmap=plt.cm.gray_r, vmin=lims[0], vmax=lims[1])
        ax.set_xlabel('Trace')
        ax.set_ylabel('Return number')

    if out_fn is not None:
        plt.savefig(out_fn)
    return fig, ax


def filter_data(data, h_hpass=0.01, h_lpass=0.1, v_hpass=0.1, v_lpass=0.2, order=2):
    """Butterworth filter int he high and/or low directions

    Parameters
    ----------
    h_hpass, h_lpass: high and low pass values for the horizontal filter (if None don't filter)
    v_hpass, v_lpass: high and low pass values for the vertical filter (if None don't filter)


    Returns
    -------
    dat: numpy.ndarray
        An array with the same size data containing the filtered data

    """
    dat = data.copy()

    # horizontal high and low passes
    if h_hpass is not None:
        b, a = signal.butter(order, h_hpass, btype='high')
        dat = signal.filtfilt(b, a, dat, axis=1)
    if h_lpass is not None:
        b, a = signal.butter(order, h_lpass, btype='low')
        dat = signal.filtfilt(b, a, dat, axis=1)

    # vertical high and low passes
    if v_lpass is not None:
        b_hor, a_hor = signal.butter(order, v_lpass)
        dat = signal.filtfilt(b_hor, a_hor, dat, axis=0)
    if v_hpass is not None:
        b_hor, a_hor = signal.butter(order, v_hpass, btype='high')
        dat = signal.filtfilt(b_hor, a_hor, dat, axis=0)
    return dat


def load_diel_file(fn):
    """Load dielectric constant from file

    Parameters
    ----------
    fn: str
        Should be a two column csv with the depths in the left column and the values in the right
    """
    with open(fn) as fin:
        lines = fin.readlines()

    # pad this with some values so we don't go out of bounds
    lines = list(map(lambda x: list(map(float, x.replace('\ufeff', '').rstrip('\n\r').split(','))), lines))
    diels = np.array([[0.0, lines[0][1]]] + lines + [[np.inf, lines[-1][1]]])
    return diels


def apply_diel(diel, data, layer_time_increment=None):
    """Apply some kind of dielectric model to data

    Parameters
    ----------
    diel: float, numpy.ndarray, str, or None
        diel can be a 2-column csv with depths and dielectric constants, an array with depths and dielectric constants, a number giving a depth-constant dielectric constant, or None, which leaves things in two-way travel time. If layer_time_increment is None, we always stay in TWTT
    data: numpy.ndarray
        The data that we are rescaling
    layer_time_increment: float
        size of a vertical increment in seconds

    Returns
    -------
    y: nx1 array
        The depths of each vertical pixel
    diels: mx2 array
        Depth and dielectric constant pairs
    """
        
    if (type(diel) == str or type(diel) == np.array) and (layer_time_increment is not None):
        # two depth-variable dielectric constants
        if type(diel) == str:
            diels = load_diel_file(diel)
        else:
            if diel.shape[1] == 2:
                diels = diel
            else:
                raise ValueError('Second dimension of dielectric array must be two')
        y = np.zeros((data.shape[0], ))

        # We now need to iterate layer by layer
        current_layer = 0
        current_dist = 0.
        y[0] = current_dist
        for i in range(1, len(y)):
            remaining_inc = layer_time_increment
            # this is complicated in case we get a sparsely sampled file where layers get skipped
            while current_layer < diels.shape[0]:
                if current_dist + remaining_inc * 3.0e8 / 2. / np.sqrt(diels[current_layer, 1]) <= diels[current_layer, 0]:
                    y[i] = current_dist + remaining_inc * 3.0e8 / 2. / np.sqrt(diels[current_layer, 1])
                    current_dist = y[i]
                    break
                else:
                    dist_this_layer = diels[current_layer, 0] - current_dist
                    time_for_this_layer = dist_this_layer / 3.0e8 * 2. * np.sqrt(diels[current_layer, 1])
                    remaining_inc -= time_for_this_layer
                    current_dist = diels[current_layer, 0]
                    current_layer += 1
        y = np.flipud(y)
    elif (diel is not None) and (layer_time_increment is not None):
        # Depth-constant diel
        incremental_travel_time = layer_time_increment * 3.0e8 / np.sqrt(diel) / 2.
        y = np.flipud(np.arange(data.shape[0])) * incremental_travel_time
        diels = np.array([[0., diel], [len(y), diel]])
    else:
        # stay in twtt
        y = np.flipud(np.arange(data.shape[0]))
        diels = np.array([[0., 1], [len(y), 1]])
    return y, diels


def zero_surface(data, method='max'):
    """Get rid of returns above the surface

    Parameters
    ----------
    data: nxm numpy.ndarray
        The radar data array
    method: str, optional
        How we detect the surface. Currently, only by the maximum amplitude return

    Returns
    -------
    data: (n-ind)xm numpy.ndarray
        The subset of the data below the surface
    ind: int
        The index of the surface pixel used as the cutoff
    """

    if method == 'max':
        avg_val = np.nanmean(data, axis=1)
        max_ind = np.argmax(avg_val)
    else:
        raise ValueError('Can only do max now')
    return data[max_ind:, :], max_ind
