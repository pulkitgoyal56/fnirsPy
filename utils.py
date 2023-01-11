#!/usr/bin/env python
# coding: utf-8

"""
Miscellaneous functions and class for fNIRS data processings (mostly in context of using MNE).
"""

import re

def is_short_channel(channel):
    """_summary_

    Parameters
    ----------
    channel : str
        Channel name

    Returns
    -------
    bool
        True if channel is short.
    """
    return re.compile(r'S(\d+)_D\1').match(channel)

def find_short_channels(channels):
    """Find short channels.

    Parameters
    ----------
    channels : array-like
        List of channels.

    Returns
    -------
    list
        List of short channels.
    list
        List of indices of short channels.
    """
    sc_dict = {i: ch for i, ch in enumerate(channels) if is_short_channel(ch)}
    return list(sc_dict.values()), list(sc_dict.keys())

def find_long_channels(channels):
    """Find long channels.

    Parameters
    ----------
    channels : array-like
        List of channels.

    Returns
    -------
    list
        List of long channels.
    list
        List of indices of long channels.
    """
    lc_dict = {i: ch for i, ch in enumerate(channels) if not is_short_channel(ch)}
    return list(lc_dict.values()), list(lc_dict.keys())

def select_best_wavelengths(wavelengths, *args):
    pass

if __name__ == '__main__':
    pass