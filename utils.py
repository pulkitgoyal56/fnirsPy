#!/usr/bin/env python
# coding: utf-8

"""
Miscellaneous functions for fNIRS data processing (mostly in context of using MNE).
"""

# File Path Manipulation
import pathlib

# Additional Iteration Utilities
from itertools import compress

# Additional Function Utilities
from functools import partial

# Regex
import re

# Array
import numpy as np

# Additional Configuration/Metadata
import constants

# MNE Constants
from mne.defaults import HEAD_SIZE_DEFAULT


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
    # TODO: Automate wavelength selection based on absorption spectra of hemoglobin.
    pass

def has_location(source, pos):
    match source:
        case pathlib.PosixPath() | str(): # source is a filename
            with open(source) as file:
                for words in file:
                    if pos == re.split('\t| +|\n', words)[0]:
                        return True
                return False
        case _: # source is a montage
            if pos in ('nasion', 'lpa', 'rpa'):
                return source.get_positions()[pos] is not None
            else:
                return pos in source.get_positions()['ch_pos']


def get_location(source, pos):
    match source:
        case pathlib.PosixPath() | str(): # source is a filename
            def _get_line_number(word, file):
                file.seek(0)
                for i, words in enumerate(file, 1):
                    if word == re.split('\t| +|\n', words)[0]:
                        return i
            def _get_word(num, file):
                file.seek(0)
                for i, words in enumerate(file, 1):
                    if num == i:
                        return words
            with open(source) as file:
                line = _get_line_number('Positions', file) - _get_line_number('Labels', file) + _get_line_number(pos, file)
                return list(map(float, re.split('\t| +', _get_word(line, file).rsplit('\n')[0])))
        case _: # source is a montage
            if pos in ('nasion', 'lpa', 'rpa'):
                return source.get_positions()[pos]
            else:
                return source.get_positions()['ch_pos'][pos]

def get_transformation(montage, reference=constants.DEFAULT_REFERENCE_LOCATIONS, scale=1/HEAD_SIZE_DEFAULT):
    """Transform montage based on expected location of reference points."""
    available_pos = list(compress(reference, map(partial(has_location, montage), reference)))

    if len(available_pos) < 4:
        raise ValueError(f'At least 4 points are requred to get complete transformation.')

    target = np.array([loc for pos, loc in reference.items() if pos in available_pos])
    target = np.c_[target, np.ones((len(available_pos), 1))]

    base = np.array(list(map(partial(get_location, montage), available_pos))) * scale
    base = np.c_[base, np.ones((len(available_pos), 1))]

    trans = np.linalg.pinv(base) @ target

    return trans.T

if __name__ == '__main__':
    pass
