#!/usr/bin/env python
# coding: utf-8

"""
NIRS class for fNIRS data processing, using MNE, for custom data.
"""

# File R/W
import os
# File Path Manipulation
import pathlib

# Regex
import re

# Array
import numpy as np
# Table
import pandas as pd

# Plotting
import matplotlib.pyplot as plt

# Advanced Computations
import scipy as sc

# Statistics
import statsmodels.api as sm

# Neurophysiological Data Analysis
import mne
import mne_nirs # esp. for fNIRS

# Misc. (custom)
import utils

# Read Data Specific Configuration/Metadata Files
import tomli

# Additional Configuration/Metadata
import constants

# Custom Functions
import utils 

class NIRS:
    def __init__(self, data_dir=constants.DATA_DIR, project=constants.PROJECT, device=constants.DEVICE) -> None:
        self.DATA_DIR = data_dir
        self.PROJECT = project
        self.DEVICE = device

    def backlight_removal(self):
        """Backlight removal."""
        raw_data = self.raw.get_data()
        time = self.raw.times

        regressors = sm.tools.tools.add_constant(np.c_[time, time**2, time**3]) # Timestamp (^1, ^2, ^3)
        
        fitted_backlight = np.apply_along_axis(
            lambda raw_backlight_ch: sm.RLM(raw_backlight_ch, regressors).fit().fittedvalues, 
            1, self.raw_backlight
        )
        
        # Subtract predicted backlight signal from raw data of all wavelengths
        corrected_data = raw_data - np.repeat(fitted_backlight, int(len(raw_data)/len(fitted_backlight)), axis=0)

        # DataFrame with backlight-removed signal intensities (the rows are chunked in groups of `N_CHANNELS` (number of used channels),
        # i.e. first 9 rows = recoding 1; second 9 rows = recording 2; ...)
        self.raw = mne.io.RawArray(corrected_data, self.raw.info)
        
        return self.raw

    def read_config(self, config_file_path):
        """Read additional configuration."""
        self.config_file_path = pathlib.Path(config_file_path).with_suffix('.toml')

        with open(self.config_file_path, 'rb') as f:
            self.config = tomli.load(f)
            
            self.S_D = utils.hex_to_dec(self.config['S_D'])
            self.WAVELENGTHS_PICKED = self.config['WAVELENGTHS_PICKED']
            self.CH_UNUSED = set(self.config['CH_UNUSED'])

    def read_raw_fif(self, raw_file_path, config_file_path=None, remove_backlight=True):
        self.raw_file_path = raw_file_path.parent / pathlib.Path(raw_file_path.stem.split('.')[0]).with_suffix('.raw.fif')

        if config_file_path:
            self.read_config(config_file_path)

        # Read '.raw.fif' file
        raw_fif = mne.io.read_raw_fif(raw_file_path, preload=True)
        # raw_fif.crop(tmin=120) # Delete first 60s for this dataset (idle data)
        
        # Wavelengths available
        self.WAVELENGTHS = pd.unique([int(ch.split()[1]) for ch in raw_fif.ch_names])
        
        # Total number of wavelengths
        self.N_WAVELENGTHS_T = len(self.WAVELENGTHS)
        
        # Number of wavelengths used in analysis
        self.N_WAVELENGTHS = len(self.WAVELENGTHS_PICKED) # 2
        
        # # Maximum number of channels
        # self.M_CHANNELS = int(len(raw_fif.ch_names) / self.N_WAVELENGTHS_T)  # per wavelength # self.config['N_PROBES'] ** 2
        
        # # Number of probes
        # self.M_PROBES = int(np.sqrt(self.M_CHANNELS) / self.self.N_HEMISPHERES) # per hemisphere # self.config['N_PROBES']

        # Source-Detector pairs (all)
        self.S_D = pd.unique([ch.split()[0] for ch in raw_fif.ch_names])

        # Short Channels
        self.CH_SHORT = set(utils.find_short_channels(self.S_D)[1]) - self.CH_UNUSED
        
        # Long Channels
        self.CH_LONG = set(utils.find_long_channels(self.S_D)[1]) - self.CH_UNUSED
        
        # Used Channels
        self.CH_USED = self.CH_SHORT.union(self.CH_LONG)

        # Number of channels
        self.N_CHANNELS = len(self.CH_USED) # per wavelength # == `int(len(raw_fif.ch_names) / N_WAVELENGTHS_T)`

        # Drop unused channels
        raw_fif.drop_channels(raw_fif.ch_names[list(self.CH_UNUSED)])

        # Names of the wavelength specific channels (only used channels)
        self.CH_NAMES = raw_fif.ch_names

        # Add custom wavelength metadata to MNE Raw Object (as private property)
        self.WAVELENGTHS_PICKED = self.WAVELENGTHS_PICKED

        # Read '-backlight.raw.fif' file
        # Backlight intensities (for used channels only)
        if remove_backlight:
            backlight_file_path = raw_file_path.parent / pathlib.Path(raw_file_path.stem.split('.')[0] + '-backlight').with_suffix('.raw.fif')
            self.raw_backlight = mne.io.read_raw_fif(backlight_file_path, preload=True).get_data()

        self.raw = raw_fif
        
        return self.raw

    def read_raw_csv(self, raw_file_path, config_file_path, remove_backlight=True):
        self.raw_file_path = pathlib.Path(raw_file_path).with_suffix('.csv')
        
        if config_file_path:
            self.read_config(config_file_path)

        # Read CSV data as Pandas DataFrame
        data_pd = pd.read_csv(raw_file_path)

        # Wavelengths available (automatic extraction)
        self.WAVELENGTHS = [int(match.groups()[0]) for column in data_pd.columns if (match := re.compile(r'(\d+)\[nm\]').match(column))] # [855, 770, 810, 885]

        # Total number of wavelengths
        self.N_WAVELENGTHS_T = len(self.WAVELENGTHS)

        # Number of wavelengths used in analysis
        self.N_WAVELENGTHS = len(self.WAVELENGTHS_PICKED) # 2

        # # Number of probes
        # self.M_PROBES = self.config['N_PROBES'] # per hemisphere

        # # Maximum number of channels
        # self.M_CHANNELS = (self.N_PROBES)**2 * self.N_HEMISPHERES # per wavelength

        # Dictionary (map) of channel number to channel name
        self.CH_MAP = dict(zip(map(int, data_pd['Channel'].unique()), self.S_D))

        # Short Channels
        self.CH_SHORT = set(utils.find_short_channels(self.S_D)[1]) - self.CH_UNUSED

        # Long Channels
        self.CH_LONG = set(utils.find_long_channels(self.S_D)[1]) - self.CH_UNUSED

        # Used Channels
        self.CH_USED = self.CH_SHORT.union(self.CH_LONG)

        # Number of channels
        self.N_CHANNELS = len(self.CH_USED) # per wavelength

        # Names of the wavelength specific channels
        # N_CHANNELS * N_WAVELENGTHS_T
        self.CH_NAMES = [f'{self.CH_MAP[ch]} {wavelength}'
                    for ch in list(self.CH_MAP.keys()) if ch not in self.CH_UNUSED
                    for wavelength in self.WAVELENGTHS]

        # Remove unused channels and create a new DataFrame
        data_pd = data_pd.loc[data_pd['Channel'].isin(self.CH_USED)]

        # # Check if the number of channels are as expected (they must not be more than `N_CHANNELS`)
        # if len(data_pd['Channel'].unique()) > M_CHANNELS:
        #     raise ValueError(f'''Duplicate channels. Expected (max.) - {M_CHANNELS}; Received - {len(data_pd["Channel"].unique())}.\n
        #                          Please pick one of the duplicates and mark the others in `CH_UNUSED` in the related config file.''')

        ## Create mne.Info object
        # Type of channels -- Raw fNIRS Continuous Wave Amplitude
        # https://mne.tools/stable/glossary.html#term-data-channels
        CH_TYPES = 'fnirs_cw_amplitude' # [] * N_CHANNELS * N_WAVELENGTHS_T * N_HEMISPHERES

        # Sampling frequency (based on difference between timestamps in consecutive readings ~54ms)
        F_S = 1000/np.ptp(data_pd['Time[ms]'])*(len(data_pd) - len(self.S_D))/len(self.S_D) # 1000/54.319 ~ 18.41 Hz

        # Create MNE.Info Object
        config_csv = mne.create_info(ch_names=self.CH_NAMES, sfreq=F_S, ch_types=CH_TYPES)

        # `Manually update info object parameters for location`
        # > Manual modification is not recommended, but there doesn't seem to any other option as there are no inbuilt functions for this.  
        # > https://github.com/mne-tools/mne-python/blob/main/mne/io/meas_config.py#L2425  
        # >> __Info__: `mne.Raw.info['chs'][x]['loc']` is an array of channel 'location' of length 12.  
        # >> From investigation, it is apparent that,  
        # >>> - [0:3] is the midpoint (channel) location (= <source_location + detector_location>/2)  
        # >>> - [3:6] is the source location  
        # >>> - [6:9] is the detector location  
        # >>> - [9] is the frequency  
        # >>> - [10] seems to be always `nan`; function unknown  
        # >>> - [11] is the separation of the channel, i.e. short (0.007) or long (0.03), in m  
        for chs in config_csv['chs']:
            # For fNIRS, the 10th element corresponds to the wavelength
            # https://github.com/mne-tools/mne-python/blob/main/mne/preprocessing/nirs/nirs.py#L150
            chs['loc'][9] = float(chs['ch_name'].split()[1])

        for chs in config_csv['chs']:
            # For fNIRS, the 12th element is the channel separation, i.e. short (0.007) or long (0.03)
            # > No specific reference to this 11th index found in the MNE-Python source code
            # >> Only references to range of values ('[:]' or '[3:]') in device-spcific functions with no apparent applicability to the context here
            chs['loc'][11] = 0.007 if utils.is_short_channel(chs['ch_name']) else 0.03

        # Copy other meta data from sample *fif* file
        config_csv['device_info'] = self.DEVICE.INFO # {'type': 'fNIRS-CW', 'model': 'optoHIVE'}
        config_csv['experimenter'] = self.DEVICE.EXPERIMENTER # 'optoHIVE Team'
        # meas_date # datetime.datetime(2022, 12, 16, 14, 36, 20, 620708, tzconfig=datetime.timezone.utc)
        # file_id (== meas_id)
        # meas_id (== file_id)

        # Create Numpy Array from the corrected DataFrame and reshape it to have rows corresponding to time-Warying signal for all channel and picked wavelength combinations
        # 'CH_USED x N_WAVELENGTHS_T' rows; each corresponding in order to `CH_NAMES`
        data_np = np.array(data_pd[data_pd.columns[-self.N_WAVELENGTHS_T:]]).reshape(-1, self.N_CHANNELS * self.N_WAVELENGTHS_T).T

        ## Create MNE.IO.Raw Object
        raw_csv = mne.io.RawArray(data_np, config_csv)

        # Add custom wavelength metadata to MNE Raw Object (as private property)
        self.WAVELENGTHS_PICKED = self.WAVELENGTHS_PICKED

        # Backlight intensities (for used channels only)
        if remove_backlight:
            self.raw_backlight = data_pd['BL'].to_numpy().reshape(-1, self.N_CHANNELS).T 

        self.raw = raw_csv

        return self.raw

    def read_raw(self, raw_file_path, config_file_path=None, remove_backlight=True):
        raw_file_path = pathlib.Path(raw_file_path)

        if raw_file_path.suffix == '':
            # Check if file with file name exists, with any extention
            if (raw_file_path := list(raw_file_path.parent.glob(f'{raw_file_path.stem}.*'))):
                raw_file_path = raw_file_path[0]
            else:
                raise FileNotFoundError(f'No file of the sort `{raw_file_path}.*`')

        match raw_file_path.suffix:
            case '.csv':
                return self.read_raw_csv(raw_file_path, config_file_path, remove_backlight)
            case '.fif':
                return self.read_raw_fif(raw_file_path, config_file_path, remove_backlight)
            case other:
                raise ValueError(f'Unsupported fNIRS file format - {other}')

    def read_annotation(self, annotation_file_path):
        self.annotation_file_path = pathlib.Path(annotation_file_path).with_suffix('.mat')
        pass

    def read_montage(self, montage_file_path):
        self.montage_file_path = pathlib.Path(montage_file_path).with_suffix('.elc')
        pass

    def read(self, subject_id, session, run):
        base_dir = pathlib.Path(self.DATA_DIR, self.PROJECT, f'sub-{subject_id}', f'ses-{session}')
        
        raw_file_path = base_dir / (f'sub-{subject_id}_ses-{session}_run-{run}_fnirs')
        annotation_file_path = base_dir / (f'sub-{subject_id}_ses-{session}_run-{run}_events.mat')
        montage_file_path = base_dir / (f'sub-{subject_id}_ses-{session}_optodes.elc')
        config_file_path = base_dir / (f'sub-{subject_id}_ses-{session}_config.toml')

        self.read_config(config_file_path)
        self.read_raw(raw_file_path)
        self.read_annotation(annotation_file_path)
        self.read_montage(montage_file_path) 
    
    def process(self, funcs):
        if type(funcs) == list:
            for func in funcs:
                self.process(func)
        else:
            self.raw = func(self)

    @staticmethod    
    def save(savepoints):
        """Saves objects passed in a dictionary."""
        def wrapper(label):
            def subwrapper(raw):
                savepoints[label] = raw
            return subwrapper
        return wrapper

    @staticmethod    
    def wrap(func, *args, **kwargs):
        def wapper(self):
            return func(self.raw, *args, **kwargs)

if __name__ == '__main__':
    pass
