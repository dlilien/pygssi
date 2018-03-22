#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dal22@uw.edu>
#
# Distributed under terms of the MIT license.

"""
A bunch of radar-related conversions
"""
import numpy as np
from datetime import datetime


def to_date(bin, le=True):
    """Convert the GSSI date format to a :class:`datetime.datetime` object

    Parameters
    ----------
    bin: bits
        The gssi binary date

    Returns
    -------
    time: :class:`datetime.datetime`
        The time from the header

    """
    def bits(bytes):
        for b in bytes:
            for i in range(8):
                yield (int(b) >> i) & 1

    def bit_to_int(bits):
        return sum([(2 ** i) * bit for i, bit in enumerate(bits)])

    bit = [b for b in bits(bin)]
    return datetime(bit_to_int(bit[25:32]) + 1980, bit_to_int(bit[21:25]), bit_to_int(bit[16:21]), hour=bit_to_int(bit[11:16]), minute=bit_to_int(bit[5:11]), second=bit_to_int(bit[0:5]) * 2)


def tt_to_m_variable_arr(diel_array, tt_arr, conv_to_sec=1.0e-9):
    """Convert two-way travel time to distance using a variable dielectric constant

    Note that this function is very much not elegant, but I think that it should not be too slow unless you dielectric array is very large

    Parameters
    ----------
    diel_array: numpy.ndarray
        An nx2 array where the first column is the depth and the second column is the dielectric constant
    tt_array: numpy.ndarray
        An mxn array of two-way travel times (typically mx1)
    conv_to_sec: float, optional
        Conversion factor to get from travel-time to seconds (default is assume in nanoseconds, which should be GSSI standard)
    """
    y = np.zeros(tt_arr.shape)
    current_layer = 0
    current_dist = 0.

    # this is complicated in case we get a sparsely sampled file where layers get skipped
    while current_layer < diel_array.shape[0]:
        bool_arr = np.logical_and(y == 0, current_dist + tt_arr * conv_to_sec * 3.0e8 / 2. / np.sqrt(diel_array[current_layer, 1]) <= diel_array[current_layer, 0])
        y[bool_arr] = (current_dist + tt_arr * conv_to_sec * 3.0e8 / 2. / np.sqrt(diel_array[current_layer, 1]))[bool_arr]

        dist_this_layer = diel_array[current_layer, 0] - current_dist
        time_for_this_layer = dist_this_layer / conv_to_sec / 3.0e8 * 2. * np.sqrt(diel_array[current_layer, 1])
        tt_arr -= time_for_this_layer
        current_dist = diel_array[current_layer, 0]
        current_layer += 1
    return y


def gained_decibels(gainpoints, data, p1=None, recenter=True):
    """Apply gain and convert to decibels

    Parameters
    ----------
    gainpoints: str
        A string of comma-separated floats for the gainpoints. These get evenly distributed across depths. They are applied exponentially
    data: numpy.ndarray
        The sample data
    p1: float, optional
        The "zero" value for the scale. If none, use the mean of `data`
    recenter: bool, optional
        If true (default), subtract the mean after gaining to center the data on zero
    """
    def apply_gain(gainpoints, data):
        gp = np.atleast_2d(np.array(list(map(float, gainpoints.split(','))))).transpose()
        space = data.shape[0] // (len(gp) - 1)
        pad_size = data.shape[0] - (len(gp) - 1) * space
        gain = np.atleast_2d(np.hstack([np.linspace(gp[i], gp[i + 1], num=space) for i in range(len(gp) - 1)] + [gp[-1] for i in range(pad_size)])).transpose()
        # Now we convert the gain to decibels
        gain = 10. ** (gain / 10.)
        data *= gain
        return data
    data2 = data_to_db(data, p1=p1)
    data2 = apply_gain(gainpoints, data2)
    if recenter:
        data2 -= data2[data2.shape[0] // 2:, :].mean()
    return data2


def data_to_db(data, p1=None):
    """Convert data to decibels

    Parameters
    ----------
    data: numpy.ndarray
        The data you want to convert
    p1: float, optional
        The "zero" value for the scale. If none, use the mean of `data`
    """
    if p1 is None:
        shp = data.shape
        p1 = data[shp[0] // 2:, :].mean()
    return 10. * np.log(data / p1)
