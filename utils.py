#!/usr/bin/env python
# coding: utf-8

"""
Miscellaneous functions for fNIRS data processing (mostly in context of using MNE).
"""

# Additional Inbuilt Utilities
from itertools import compress
# Regex
import re

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

is_long_channel = lambda channel: not is_short_channel(channel)

def find_short_channels(channels):
    """Find short channels from names.

    Parameters
    ----------
    channels : array-like
        List of channel names.

    Returns
    -------
    list
        List of short channels (names).
    list
        List of indices of short channels.
    """
    return list(filter(is_short_channel, channels)), list(compress(range(len(channels)), map(is_short_channel, channels)))

def find_long_channels(channels):
    """Find long channels from names.

    Parameters
    ----------
    channels : array-like
        List of channel names.

    Returns
    -------
    list
        List of long channels (names).
    list
        List of indices of long channels.
    """
    return list(filter(is_long_channel, channels)), list(compress(range(len(channels)), map(is_long_channel, channels)))

def select_best_wavelengths(wavelengths, *args):
    pass

if __name__ == '__main__':
    pass
