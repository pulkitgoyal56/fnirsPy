#!/usr/bin/env python
# coding: utf-8

"""
Miscellaneous functions and classes for fNIRS data processings (mostly in context of using MNE).
"""

import re
from itertools import compress

def dec_to_hex(channels):
    """Converts channel IDs from decimal to hexadecimal.

    Parameters
    ----------
    channels : array-like
        List of channel names with decimal IDs.

    Returns
    -------
    array-like
        List of channel names with hexadecimal IDs.
    """
    def convert(channel):
        s, d, r = re.compile(r'S(\d+)_D(\d+)(.*)').match(channel).groups()
        return f'S{int(s):X}_D{int(d):X}{r}'
    
    return list(map(convert, channels))

def hex_to_dec(channels):
    """Converts channel IDs from hexadecimal to decimal.

    Parameters
    ----------
    channels : array-like
        List of channel names with hexadecimal IDs.

    Returns
    -------
    array-like
        List of channel names with decimal IDs.
    """
    def convert(channel):
        s, d, r = re.compile(r'S([0-9A-F]+)_D([0-9A-F]+)(.*)').match(channel).groups()
        return f'S{int(s, base=16)}_D{int(d, base=16)}{r}'
    
    return list(map(convert, channels))

def is_short_channel(channel):
    """Check if channel is short based on its name.

    Parameters
    ----------
    channel : str
        Channel name.

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