#!/usr/bin/env python
# coding: utf-8

"""
Miscellaneous functions for fNIRS data processing (mostly in context of using MNE).
"""

# File Path Manipulation
import pathlib

# Additional Iteration Utilities
from itertools import compress
from collections import Counter

# Additional Function Utilities
from functools import partial

# Regex
# TIP - Use regexr.com to create/test/learn Regex
import re

# Array
import numpy as np

# Additional Configuration/Metadata
import constants

# MNE Constants
from mne.defaults import HEAD_SIZE_DEFAULT


def dec_to_hex(ch_names):
    """Converts channel IDs from decimal to hexadecimal.

    Parameters
    ----------
    ch_names : array-like
        List of channel names with decimal IDs.

    Returns
    -------
    list
        List of channel names with hexadecimal IDs.
    """
    def convert(ch_name):
        s, d, r = re.compile(r'S(\d+)_D(\d+)(.*)').match(ch_name).groups()
        return f'S{int(s):X}_D{int(d):X}{r}'

    return list(map(convert, ch_names))

def hex_to_dec(ch_names):
    """Converts channel IDs from hexadecimal to decimal.

    Parameters
    ----------
    ch_names : array-like
        List of channel names with hexadecimal IDs.

    Returns
    -------
    list
        List of channel names with decimal IDs.
    """
    def convert(ch_name):
        s, d, r = re.compile(r'S([0-9A-F]+)_D([0-9A-F]+)(.*)').match(ch_name).groups()
        return f'S{int(s, base=16)}_D{int(d, base=16)}{r}'

    return list(map(convert, ch_names))

def get_s_d(ch_names):
    """Gets the unique source detector names without wavelength/chromophore labels.

    Parameters
    ----------
    ch_names : array-like
        List of channel names with wavelength/chromophore labels.

    Returns
    -------
    list
        Ordered list of unique channel names without wavelength/chromophore labels.
    """
    return list(dict.fromkeys(map(lambda ch_name: ch_name.split()[0], ch_names)))

def is_short_channel(ch_name):
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
    return re.compile(r'S(\d+)_D\1').match(ch_name)

is_long_channel = lambda ch_name: not is_short_channel(ch_name)

def is_channel_type(type, ch_name):
    """Check if channel is of given type based on its name.

    Parameters
    ----------
    type: str
        Channel type.
    channel : str
        Channel name.

    Returns
    -------
    bool
        True if channel is of given type.
    """
    return re.compile(fr'S\d+_D\d+ {type}').match(ch_name)

def find_short_channels(ch_names):
    """Find short channels from names.

    Parameters
    ----------
    ch_names : array-like
        List of channel names.

    Returns
    -------
    list
        List of short channels (names).
    list
        List of indices of short channels.
    """
    return list(filter(is_short_channel, ch_names)), list(compress(range(len(ch_names)), map(is_short_channel, ch_names)))

def find_long_channels(ch_names):
    """Find long channels from names.

    Parameters
    ----------
    ch_names : array-like
        List of channel names.

    Returns
    -------
    list
        List of long channels (names).
    list
        List of indices of long channels.
    """
    return list(filter(is_long_channel, ch_names)), list(compress(range(len(ch_names)), map(is_long_channel, ch_names)))

def find_channels_type(type, ch_names):
    """Find channels of given type from names.

    Parameters
    ----------
    type: str
        Channel type.
    ch_names : array-like
        List of channel names.

    Returns
    -------
    list
        List of channels (names) of given type.
    list
        List of indices of these channels.
    """
    return list(filter(partial(is_channel_type, type), ch_names)), list(compress(range(len(ch_names)), map(partial(is_channel_type, type), ch_names)))

def select_best_wavelengths(wavelengths, *args):
    # TODO: Automate wavelength selection based on absorption spectra of hemoglobin.
    pass

def find_ch_pairs(ch_names, channels):
    """Find the other channels with the same source-detector pair as the queried channels.

    Parameters
    ----------
    ch_names : array-like
        List of channel names.
    channels : array-like, str or int
        List of channel names or ids to be paired.

    Returns
    -------
    list
        Names/ids of remaining channels with same source-detector pair.
    """
    match next(iter(channels)):
        case int():
            id = True
            channels = [ch_names[ch] for ch in channels]
        case str():
            id = False

    return [ch if id else ch_name for ch, ch_name in enumerate(ch_names) if ch_name not in channels and get_s_d([ch_name])[0] in get_s_d(channels)]

def find_ch_unpaired(channels, ch_names=None):
    """Find the channels that are not paired.

    Parameters
    ----------
    channels : array-like, str or int
        List of channel names or ids to be filtered.
    ch_names : array-like, optional
        List of channel names.
        Required if channels is a list of ids.

    Returns
    -------
    list
        Names/ids of channels with same source-detector pair.
    """
    match next(iter(channels)):
        case int():
            id = True
            channels = [ch_names[ch] for ch in channels]
        case str():
            id = False

    return [ch_names.index(ch_name) if id else ch_name for ch_name in channels if Counter(map(lambda ch_name: ch_name.split()[0], channels))[get_s_d([ch_name])[0]] == 1]

def has_location(source, pos):
    match source:
        case pathlib.PurePath() | str(): # source is a filename
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
        case pathlib.PurePath() | str(): # source is a filename
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
