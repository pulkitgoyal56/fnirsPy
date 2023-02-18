#!/usr/bin/env python
# coding: utf-8

"""
NIRS class for fNIRS data processing, using MNE, for custom data.
"""

# Annotations for type checking
from __future__ import annotations

# Logging and warnings
import logging
import warnings

# File R/W
import os

# File Path Manipulation
import pathlib

# Additional Inbuilt Utilities
import itertools
import functools

# Regex
import re

# Array
import numpy as np
# Table
import pandas as pd

# MATLAB v7.3 File Format
import mat73

# Advanced Computations
import scipy as sc

# Statistics
import statsmodels.api as sm

# Neurophysiological Data Analysis
import mne
import mne_nirs # esp. for fNIRS

# Neuroimaging Statistical Tools
import nilearn
import nilearn.plotting

# Plotting
import matplotlib.pyplot as plt

# Read Data Specific Configuration/Metadata Files
import tomli

# Additional Configuration/Metadata
import constants

# Custom Functions
import utils

# Custom MBLL
import mbll


class NIRS:
    def __init__(self, data_dir=constants.DATA_DIR, project=constants.PROJECT, device=constants.DEVICE) -> None:
        """Initialize NIRS object."""
        self._DATA_DIR = data_dir
        self._PROJECT = project
        self._DEVICE = device

        self._TIME_DRIFT_FACTOR = 1.0

    def __len__(self):
        """Number of channels (excl. bad channels)."""
        return len(self.raw.ch_names) - len(self.raw.info['bads'])

    @property
    def shape(self):
        """Shape of data (excl. bad channels)."""
        n_chs, n_t = self.raw.get_data().shape
        return n_chs - len(self.raw.info['bads']), n_t

    def __attr(self, attribute, value):
        """Modify the object attribute if value is given, otherwise return the object attribute value"""
        if value is not None:
            self.__setattr__(attribute, value)
        return self.__getattribute__(attribute)

    def pick_wavelengths(self, wavelengths_picked=None, **kwargs):
        """Pick wavelengths."""
        wavelengths_picked = self.__attr('WAVELENGTHS_PICKED', wavelengths_picked)
        self.wavelengths = wavelengths_picked

        # Indices of all the channels available (beware, these are not the same as the initial channel numbers!)
        # picks = mne.pick_types(self.raw.info, meg=False, fnirs=True) # Select channels with picked wavelengths

        # Pick long channels (for picked wavelength)
        # self.raw.pick([ch_name for ch_name in utils.find_long_channels(self.raw.ch_names)[0] if int(ch_name.split()[1]) in self.WAVELENGTHS_PICKED])

        self.raw.pick([ch_name for ch_name in self.raw.ch_names if int(ch_name.split()[1]) in self.wavelengths])

    def set_bad(self, bad_channels, *, overwrite=False):
        """Set bad channels."""
        if isinstance(next(iter(bad_channels)), int):
            bad_channels = [self.raw.ch_names[ch] for ch in bad_channels]

        self.raw.info['bads'] = [ch_name for ch_name in self.raw.ch_names if ch_name in bad_channels or not overwrite and ch_name in self.raw.info['bads']]

    def read_config(self, config_file_path, **kwargs):
        """Read additional configuration."""
        self.config_file_path = pathlib.Path(config_file_path).with_suffix('.toml')

        with open(self.config_file_path, 'rb') as f:
            self.CONFIG = tomli.load(f)

            self.WAVELENGTHS_PICKED = self.CONFIG['WAVELENGTHS_PICKED']
            self.S_D = utils.hex_to_dec(self.CONFIG['S_D'])
            self.T_EXP_START = float(self.CONFIG['T_EXP_START'])
            self.T_EPOCH_START = float(self.CONFIG['T_EPOCH_START'])
            self.T_EPOCH_END = float(self.CONFIG['T_EPOCH_END'])
            self.T_BASELINE_START = float(self.CONFIG['T_BASELINE_START'])
            self.T_BASELINE_END = float(self.CONFIG['T_BASELINE_END'])

    @property
    def s_d(self):
        """Get source detector pairs (incl. bad channels)."""
        # # Do not count channels marked as bad, if all frequencies/chromophores for that channel are marked bad.
        # # If any of the chromophore is not marked bad, count it still!
        # return list(dict.fromkeys([ch_name.split(' ')[0] for ch_name in self.raw.ch_names if ch_name not in self.raw.info['bads']]))

        return utils.get_s_d(self.raw.ch_names)

    @property
    def good_ch_names(self):
        """Get names of good channels in order."""
        # list(set(self.raw.ch_names) - set(self.raw.info['bads']))
        return list(dict.fromkeys([ch_name for ch_name in self.raw.ch_names if ch_name not in self.raw.info['bads']]))

    @property
    def bad_ch_names(self):
        """Get names of bad channels."""
        return self.raw.info['bads']

    def correct_time(self, correction_factor=constants.DEVICE.TIME_DRIFT_FACTOR, **kwargs):
        match correction_factor:
            case 'auto':
                correction_factor = self.DUR['rec'] / self.DUR['exp']
                logging.info(f'''Using Auto Correction Factor - {correction_factor}''')
            case 'default' | True:
                correction_factor = self._DEVICE.TIME_DRIFT_FACTOR
            case False | None:
                correction_factor = 1.0

        correction_factor /= self._TIME_DRIFT_FACTOR

        self.F_S *= correction_factor
        self.T_REC_START /= correction_factor
        self.T_REC_END /= correction_factor
        self.DUR['rec'] /= correction_factor

        self._TIME_DRIFT_FACTOR *= correction_factor

        # Update sampling frequency
        # Note - To update certain attributes of the mne.Info object, the state has to manually 'unlocked'
        #      - This is not recommended, but the alternative to recreate the info object for one change is unacceptable
        #      - If there are any sync issues with info object, this might be the where to investigate but it likely won't happen
        #      - The attributes that can be updated without 'unlocking' are
        #      - > `info['bads']`, `info['description']`, `info['device_info'`], `info['dev_head_t']`, `info['experimenter']`,
        #      - > `info['helium_info']`, `info['line_freq']`, `info['temp']`, and `info['subject_info']`
        #      - "All other entries should be considered read-only, though they can be modified by various MNE-Python functions or methods
        #      - (which have safeguards to ensure all fields remain in sync)."
        #      - See - https://mne.tools/stable/generated/mne.Info.html#mne.Info
        self.raw.info._unlocked = True
        self.raw.info['sfreq'] = self.F_S
        self.raw.info._unlocked = False

    @property
    def tdf(self):
        """Getter method for `_TIME_DRIFT_FACTOR`."""
        return self._TIME_DRIFT_FACTOR

    @tdf.setter
    def tdf(self, correction_factor):
        """Setter method in case the time drift factor is manually modified."""
        self.correct_time(correction_factor)

    def read_raw_fif(self, raw_file_path, config_file_path=None, *, backlight=True, **kwargs):
        """Read .fif file and its accompanying backlight file."""
        raw_file_path = pathlib.Path(raw_file_path)
        self.raw_file_path = raw_file_path.parent / pathlib.Path(raw_file_path.stem.split('.')[0]).with_suffix('.raw.fif')

        if config_file_path:
            self.read_config(config_file_path)

        # Read '.raw.fif' file and create mne.raw object
        raw_fif = mne.io.read_raw_fif(self.raw_file_path, preload=True)

        # Wavelengths available
        self.wavelengths = self.WAVELENGTHS = pd.unique([int(ch_name.split()[1]) for ch_name in raw_fif.ch_names])

        # # Maximum number of channels
        # self.M_CHANNELS = int(len(raw_fif.ch_names) / len(wavelengths))  # per wavelength # int(self.CONFIG['N_PROBES']) ** 2

        # # Number of probes
        # self.M_PROBES = np.ceil(np.sqrt(self.M_CHANNELS) / int(self.CONFIG['N_HEMISPHERES']))  # per hemisphere # int(self.CONFIG['N_PROBES'])

        # Source-Detector Pairs (all)
        self.S_D = utils.get_s_d(raw_fif.ch_names)

        # Used Channels
        self.S_D_USED = [s_d for i, s_d in enumerate(self.S_D) if i not in self.CONFIG['S_D_UNUSED']] # per wavelength

        # Short Channels
        # S_D_SHORT = utils.find_short_channels(self.S_D_USED)[1]

        # Long Channels
        # S_D_LONG = utils.find_long_channels(self.S_D_USED)[1]

        # Names of the wavelength specific channels (only used channels)
        self.CH_NAMES = raw_fif.ch_names

        # Drop unused channels
        raw_fif.drop_channels([ch_name for ch_name in self.CH_NAMES if utils.get_s_d([ch_name])[0] not in self.S_D_USED])

        assert len(raw_fif.ch_names) == len(self.S_D_USED) * len(self.WAVELENGTHS)

        # Sampling frequency (based on difference between timestamps in consecutive readings ~54ms)
        self.F_S = raw_fif.info['sfreq']                                     # fNIRS recording frequency, in Hertz

        # Set recording start and end times
        self.T_REC_START = 0                                                 # fNIRS recording start time, in seconds
        self.T_REC_END = len(raw_fif)/self.F_S                               # fNIRS recording end time, in seconds

        # Re-set channel types in case different data is read
        self.CONFIG['CH_TYPES'] = raw_fif.info.get_channel_types()

        # Re-set other meta data
        self._DEVICE.INFO = raw_fif.info['device_info'] # {'type': 'fNIRS-CW', 'model': 'optoHIVE'}
        self._DEVICE.EXPERIMENTER = raw_fif.info['experimenter'] # 'optoHIVE Team'

        # Re-create mne.raw object
        self.raw = raw_fif

        # Read '-backlight.raw.fif' file
        # Backlight intensities (for used channels only)
        if backlight:
            backlight_file_path = raw_file_path.parent / pathlib.Path(raw_file_path.stem.split('.')[0] + '-backlight').with_suffix('.raw.fif')
            self.raw_backlight = mne.io.read_raw_fif(backlight_file_path, preload=True).get_data()

        return self.raw

    def read_raw_csv(self, raw_file_path, config_file_path=None, *, backlight=True, **kwargs):
        """Read .csv file with its accompanying backlight data."""
        self.raw_file_path = pathlib.Path(raw_file_path).with_suffix('.csv')

        if config_file_path:
            self.read_config(config_file_path)

        # Read CSV data as Pandas DataFrame
        # The rows are chunked in groups of n_s_d (number of source-detector pairs),
        # e.g. first <n_s_d> rows = recording 1; second x rows = recording 2; ...)
        data_pd = pd.read_csv(self.raw_file_path)

        # Wavelengths available (automatic extraction)
        self.wavelengths = self.WAVELENGTHS = [int(match.groups()[0]) for column in data_pd.columns if (match := re.compile(r'(\d+)\[nm\]').match(column))]

        # # Number of probes
        # self.M_PROBES = int(self.CONFIG['N_PROBES']) # per hemisphere

        # # Maximum number of channels
        # self.M_CHANNELS = (self.M_PROBES)**2 * int(self.CONFIG['N_HEMISPHERES']) # per wavelength

        # Used Channels
        self.S_D_USED = [s_d for i, s_d in enumerate(self.S_D) if i not in self.CONFIG['S_D_UNUSED']] # per wavelength

        # Short Channels
        # S_D_SHORT = utils.find_short_channels(self.S_D_USED)[1]

        # Long Channels
        # S_D_LONG = utils.find_long_channels(self.S_D_USED)[1]

        # Names of the wavelength specific channels
        self.CH_NAMES = [f'{s_d} {wavelength}' for s_d in self.S_D_USED for wavelength in self.WAVELENGTHS]

        if len(data_pd['Channel'].unique()) != len(self.S_D):
            warnings.warn(f'''The number of source-detector pairs provided ({len(self.S_D)}) is not equal to the number of channels in the data ({len(data_pd['Channel'].unique())})!''')

        # Remove unused channels and create a new DataFrame
        data_pd = data_pd.loc[data_pd['Channel'].isin([self.S_D.index(s_d) for s_d in self.S_D_USED])]

        assert len(data_pd['Channel'].unique()) == len(self.S_D_USED)

        # # Check if the number of channels are as expected (they must not be more than `M_CHANNELS`)
        # if len(data_pd['Channel'].unique()) > M_CHANNELS:
        #     raise ValueError(f'''Duplicate channels. Expected (max.) - {M_CHANNELS}; Received - {len(data_pd["Channel"].unique())}.\n
        #                          Please pick one of the duplicates and mark the others in `S_D_UNUSED` in the related config file.''')

        # Set recording start and end times
        self.T_REC_START = -data_pd['Time[ms]'].iloc[0]/1000                                # fNIRS recording start time, in seconds
        self.T_REC_END = np.ptp(data_pd['Time[ms]'])/1000                                   # fNIRS recording end time, in seconds

        # Sampling frequency (based on difference between timestamps in consecutive readings ~54ms)
        self.F_S = len(data_pd)/len(self.S_D_USED)/self.T_REC_END    # fNIRS recording frequency, in Hertz

        # Create mne.Info Object
        info_csv = mne.create_info(ch_names=self.CH_NAMES, sfreq=self.F_S, ch_types=self.CONFIG['CH_TYPES'])

        # `Manually update info object parameters for location`
        # > Manual modification is not recommended, but there doesn't seem to any other option as there are no inbuilt functions for this.
        # > https://github.com/mne-tools/mne-python/blob/main/mne/io/meas_config.py#L2425
        # >> __Info__: `mne.raw.info['chs'][x]['loc']` is an array of channel 'location' of length 12.
        # >> From investigation, it is apparent that,
        # >>> - [0:3] is the midpoint (channel) location (= <source_location + detector_location>/2)
        # >>> - [3:6] is the source location
        # >>> - [6:9] is the detector location
        # >>> - [9] is the wavelength
        # >>> - [10] seems to be always `nan`; function unknown
        # >>> - [11] is the separation of the channel, in m
        for chs in info_csv['chs']:
            # For fNIRS, the 10th element corresponds to the wavelength
            # https://github.com/mne-tools/mne-python/blob/main/mne/preprocessing/nirs/nirs.py#L150
            chs['loc'][9] = float(chs['ch_name'].split()[1])

        for chs in info_csv['chs']:
            # For fNIRS, the 12th element is the channel separation
            # > No specific reference to this 11th index found in the MNE-Python source code
            # >> Only references to range of values ('[:]' or '[3:]') in device-spcific functions with no apparent applicability to the context here
            chs['loc'][11] = self._DEVICE.SS_SEPARATION if utils.is_short_channel(chs['ch_name']) else self._DEVICE.LS_SEPARATION

        # Copy other meta data
        info_csv['device_info'] = self._DEVICE.INFO # {'type': 'fNIRS-CW', 'model': 'optoHIVE'}
        info_csv['experimenter'] = self._DEVICE.EXPERIMENTER # 'optoHIVE Team'
        # meas_date # datetime.datetime(2022, 12, 16, 14, 36, 20, 620708, tzconfig=datetime.timezone.utc)
        # file_id (== meas_id)
        # meas_id (== file_id)

        # Create Numpy Array from the corrected DataFrame and reshape it to have rows corresponding to time-warying signal for all channel and picked wavelength combinations
        # 'n_channels = n_s_d x n_wavelengths' rows; each corresponding in order to `CH_NAMES`
        data_np = np.array(data_pd[data_pd.columns[-len(self.WAVELENGTHS):]]).reshape(-1, len(self.S_D_USED) * len(self.WAVELENGTHS)).T

        # Create mne.raw object
        self.raw = mne.io.RawArray(data_np, info_csv)

        # Backlight intensities (for used channels only)
        if backlight:
            self.raw_backlight = data_pd['BL'].to_numpy().reshape(-1, len(self.S_D_USED)).T

        return self.raw

    def read_raw(self, raw_file_path, config_file_path=None, **kwargs):
        """Read raw data."""
        raw_file_path = pathlib.Path(raw_file_path)

        if raw_file_path.suffix == '':
            # Check if file with file name exists, with any extention
            if (raw_file_path := list(raw_file_path.parent.glob(f'{raw_file_path.stem}.*'))):
                raw_file_path = raw_file_path[0]
            else:
                raise FileNotFoundError(f"No file of the sort `{raw_file_path}.*`")

        match raw_file_path.suffix:
            case '.fif':
                return self.read_raw_fif(raw_file_path, config_file_path, **kwargs)
            case '.csv':
                return self.read_raw_csv(raw_file_path, config_file_path, **kwargs)
            case other:
                raise ValueError(f"Unsupported fNIRS file format - {other}")

    def read_annotation(self, annotation_file_path, **kwargs):
        """Read annotation data."""
        self.annotation_file_path = pathlib.Path(annotation_file_path).with_suffix('.mat')

        if self._PROJECT == 'Multi-object Tracking':
            # `Stages of the experiment`
            # > *\<exp\>* → **\[ *\<tri\>* = *\<wait1\>* → *\<target\>* → *\<motion\>* → *\<probe\>* → *\<feedb\>* → *\<feedbEnd\>* → {data_write()} \]** → *\<expEnd\>*
            # > *`T_REC_START`* ------ *`T_EXP_START`* == *0* ------------------------------------------------------------ *`T_EXP_END`* ------ *`T_REC_END`*

            # Load Experiment Results
            self.mat = pd.DataFrame(sc.io.loadmat(self.annotation_file_path)['blockdata'], columns=[
                'experiment_number',     # 0     # <trl.exp_number>              # experiment number
                'subject_number',        # 1     # <trl.sub_number>              # subject number
                'trial_number',          # 2     # <trl.num>                     # trial number

                'num_targets',           # 3     # <trl.numTargets(trl.num)>     # number of targets ({0, 2, 3, 4, 5})
                'probe_match',           # 4     # <trl.probeMatch(trl.num)>     # target (1) or not (0)

                'checker_pres',          # 5     # <trl.checkerPres(trl.num)>    # 0 -- checkerboard present (1) or not (0)
                'checker_side',          # 6     # <trl.checkerSide(trl.num)>    # 0 -- checkerboard display side; right (1) or left (2), or NA (0)

                'id_correct',            # 7     # <dat.IDcor>                   # correct (1), incorrect (2), other (3), none (4)

                'rt_correct',            # 8     # <dat.RTcor>                   # {∈ (trl.RTmin, trl.RTmax)} (1), {≤ trl.RTmin (2)}, {≥ trl.RTmax} (3), {> trl.RTmaxWait} (4)

                'rt',                    # 9     # <dat.RT>                      # reaction time, in seconds

                'tri_e',                 # 10    # <time.triE>                   # <tri - exp> -- trial starting time, relative to start of the experiment, in seconds
                'wait_e',                # 11    # <time.waitE>                  # <wait1 - exp> -- wait starting time, relative to the start of the experiment, in seconds
                'target_e',              # 12    # <time.targetE>                # <target - exp> -- target display time, relative to start of the experiment, in seconds
                'motion_e',              # 13    # <time.motionE>                # <motion - exp> -- motion starting time, relative to start of the experiment, in seconds
                'probe_e',               # 14    # <time.probeE>                 # <probe - exp> -- response collection time, relative to start of the experiment, in seconds

                'tri_t',                 # 15    # <time.triT>                   # <tri - tri> -- trial starting time, relative to start of the trial, in seconds
                'wait_t',                # 16    # <time.waitT>                  # <wait - tri> -- ??
                'target_t',              # 17    # <time.targetT>                # <target - tri> -- target display time, relative to start of the trial, in seconds
                'motion_t',              # 18    # <time.motionT>                # <motion - tri> -- motion starting time, relative to start of the trial, in seconds
                'probe_t',               # 19    # <time.probeT>                 # <probe - tri> -- response collection time, relative to start of the trial, in seconds

                'wait_p',                # 20    # <time.waitP>                  # <target - wait1> -- wait duration, in seconds
                'target_p',              # 21    # <time.targetP>                # <motion - target> -- target presentation duration, in seconds
                'motion_p',              # 22    # <time.motionP>                # <probe - motion> -- motion duration, in seconds
                'probe_p',               # 23    # <time.probeP>                 # <feedb - probe> -- time window for response input, in seconds
                'feedb_p'                # 24    # <time.feedbP>                 # <feedbEnd - feedb> -- response feedback report duration, in seconds
            ])

            # Create dictionary of all the durations of a trial by looking for all the columns with names ending with '_p'
            self.DUR = self.mat.filter(regex=('.*_p$')).mean().round().rename(lambda c_n: c_n[:-2]) # .astype(int)
            self.DUR['trial'] = sum(self.DUR)
            self.F_E = 1/self.DUR['trial']

            # Read experiment end time, relative to its start time, vis-à-vis its duration
            endtime_file_path = annotation_file_path.parent / pathlib.Path(annotation_file_path.stem.rsplit('_', 1)[0] + '_endtime').with_suffix('.mat')
            self.DUR['exp'] = float(sc.io.loadmat(endtime_file_path)['expEnd'])  # <expEnd - exp> -- duration of the entire experiment, in seconds

            # Set the duration of the recording
            self.DUR['rec'] = self.T_REC_END     # - self.T_REC_START            # recording duration, in seconds

            # Read experiment end time and set start time
            # self.T_EXP_START = 0               # <exp>                         # experiment start time, in seconds; offset due to trigger delay, in seconds
            self.T_EXP_END = self.T_EXP_START + self.DUR['exp']                  # experiment end time, in seconds

            # There could be time differences in the fNIRS recordings and experiment, due to fast/slow clocks of the device.
            # This can be corrected by scaling the recording times and frequencies by a correction factor.
            # The correction factor is greater than 1 if `DUR['rec']` > `DUR['exp']`, and vice versa.
            self.correct_time(**kwargs)

            # Set annotations in the Raw object
            self.raw.set_annotations(mne.Annotations(
                onset=self.T_EXP_START + self.mat['motion_e'], # - self.T_REC_START
                duration=self.mat['motion_p'],
                description=self.mat['num_targets'].astype(int) # TODO: Read alternative annotation descriptions from kwargs or introduce new `description` argument.
            ))

        elif self._PROJECT == 'Working-Memory':
            # `Stages of the experiment`
            # > *\<exp\>* → **\[ *\<tri\>* → *\<ti.ms\>* → *\<ti.mi\>* → *\<ti.mp\>* → *\<ti.ri\>* → *\<ti.iti\>* → {data_write()} \]** → *\<expEnd\>*
            # > *`T_REC_START`* ------ *`T_EXP_START`* == *0* ------------------------------------------------------------ *`T_EXP_END`* ------ *`T_REC_END`*

            # Load Experiment Results
            self.annotation_file_path = annotation_file_path.parent / pathlib.Path(annotation_file_path.stem.rsplit('_', 1)[0] + '_onsets').with_suffix('.mat')
            self.mat = pd.DataFrame(sc.io.loadmat(self.annotation_file_path)['onsets'], columns=[
                'onsets',                # 0     # <>                            # onsets
                'condition'              # 1     # <>                            # load condition
            ])

            # Create dictionary of all the durations of a trial
            self.DUR = pd.Series({
                'target': 1.0,                                                   # <ti.ms> -- memory set
                'motion': 6.0,                                                   # <ti.mi> -- maintainance interval
                'probe': 1.0,                                                    # <ti.mp> -- memory probe
                'feedb': 5.0,                                                    # <ti.ri> -- response interval
                'wait': 1.0,                                                     # <ti.ti> -- inter trial interval
            })
            self.DUR['trial'] = sum(self.DUR)
            self.F_E = 1/self.DUR['trial']

            # Read experiment end time, relative to its start time, vis-à-vis its duration
            endtime_file_path = annotation_file_path.parent / pathlib.Path(annotation_file_path.stem.rsplit('_', 1)[0] + '_endtime').with_suffix('.mat')
            self.DUR['exp'] = float(sc.io.loadmat(endtime_file_path)['expEnd'])  # <expEnd - exp> -- duration of the entire experiment, in seconds

            # Set the duration of the recording
            self.DUR['rec'] = self.T_REC_END     # - self.T_REC_START            # recording duration, in seconds

            # Read experiment end time and set start time
            # self.T_EXP_START = 0               # <exp>                         # experiment start time, in seconds; offset due to trigger delay, in seconds
            self.T_EXP_END = self.T_EXP_START + self.DUR['exp']                  # experiment end time, in seconds

            # There could be time differences in the fNIRS recordings and experiment, due to fast/slow clocks of the device.
            # This can be corrected by scaling the recording times and frequencies by a correction factor.
            # The correction factor is greater than 1 if `DUR['rec']` > `DUR['exp']`, and vice versa.
            self.correct_time(**kwargs)

            # Set annotations in the Raw object
            self.raw.set_annotations(mne.Annotations(
                onset=self.mat['onsets'], # - self.T_REC_START
                duration=self.DUR['motion'],
                description=('LOW', 'HIGH')[int(np.unique(self.mat['condition'])) - 1]  # TODO: Read alternative annotation descriptions from kwargs or introduce new `description` argument.
            ))

        # Extract events of interest
        self.events, self.event_dict = mne.events_from_annotations(self.raw)
        self.cases = list(self.event_dict)

        return self.events, self.event_dict

    def read_montage(self, montage_file_path, *, augment=True, transform=True, reference_locations=constants.DEFAULT_REFERENCE_LOCATIONS, reference=constants.DEFAULT_REFERENCE, **kwargs):
        """Read location data."""
        self.montage_file_path = pathlib.Path(montage_file_path).with_suffix('.elc')

        montage = mne.channels.read_custom_montage(self.montage_file_path, coord_frame=self.CONFIG['COORD_FRAME'], head_size=self.CONFIG['HEAD_SIZE'])
        # TIP - The montage stores location after dividing by a constant factor of order 3

        # Add missing fiducial point nasion coordinates from MNI coordinates if missing
        if augment and (montage.get_positions()['coord_frame'] == 'mri'):
            montage.add_estimated_fiducials('fsaverage', mne.datasets.sample.data_path() / 'subjects')

        if transform:
            match reference_locations:
                case 'default' | dict():
                    if reference_locations == 'default':
                        reference_locations = constants.DEFAULT_REFERENCE_LOCATIONS
                    montage.apply_trans(mne.transforms.Transform(fro=self.CONFIG['COORD_FRAME'], to=self.__attr('REFERENCE', reference),
                        trans=utils.get_transformation(montage, reference_locations,
                                scale=utils.get_location(self.montage_file_path, next(iter(montage.get_positions()['ch_pos'])))[0]/
                                        next(iter(montage.get_positions()['ch_pos'].values()))[0])))
                case str():
                    # TODO: Add other specific cases and use inbuilt transformations.
                    # see https://github.com/mne-tools/mne-python/blob/maint/1.3/mne/transforms.py#L641
                    raise ValueError(f"{reference_locations} method is not supported yet.")
                case _:
                    raise ValueError(f"Unsupported `reference_locations`.")

        # Picking wavelengths is required because MNE does not support more than two wavelengths.
        # An error will be raised when setting the montage if more than two wavelengths are used.
        self.pick_wavelengths(**kwargs)

        # montage.plot()
        self.montage = montage
        self.raw.set_montage(self.montage)

    def read(self, subject_id, session, run, **kwargs):
        """Read subject/session/run data; set annotations and set montage."""
        base_dir = pathlib.Path(self._DATA_DIR, self._PROJECT, f'sub-{subject_id}', f'ses-{session}')

        raw_file_path = base_dir / (f'sub-{subject_id}_ses-{session}_run-{run}_fnirs')
        annotation_file_path = base_dir / (f'sub-{subject_id}_ses-{session}_run-{run}_events.mat')
        montage_file_path = base_dir / (f'sub-{subject_id}_ses-{session}_optodes.elc')
        config_file_path = base_dir / (f'sub-{subject_id}_ses-{session}_config.toml')

        self.read_config(config_file_path, **kwargs)
        self.read_raw(raw_file_path, **kwargs)
        self.read_annotation(annotation_file_path, **kwargs)
        self.read_montage(montage_file_path, **kwargs)

        return self

    def process(self, *funcs):
        """Process NIRS object with function that takes the object as input and returns mne.raw instance if it's modified, else whatever."""
        for func in funcs:
            if isinstance((raw := func(self)), type(self.raw)):
                self.raw = raw

    @staticmethod
    def save(savepoints):
        """Saves objects passed in a dictionary."""
        def wrapper(label):
            def subwrapper(self):
                if isinstance(savepoints, dict):
                    savepoints[label] = self.raw.copy()
            return subwrapper
        return wrapper

    @staticmethod
    def wrap(func):
        """Wraps functions that take raw to take NIRS object."""
        @functools.wraps(func)
        def wrapper(*args, execute=True, **kwargs):
            @functools.wraps(func)
            def subwrapper(self):
                if execute:
                    return func(self.raw, *args, **kwargs)
            return subwrapper
        return wrapper

    def get_epochs(self, tmin=None, tmax=None, baseline=(None, None), reject_criteria=constants.REJECT_CRITERIA, *,
                   reject_by_annotation=True, preload=True, plot_drop_log=False, **kwargs):
        """Extract epochs."""
        tmin = self.__attr('T_EPOCH_START', tmin)
        tmax = self.__attr('T_EPOCH_END', tmax)

        baseline = (self.__attr('T_BASELINE_START', baseline[0]), self.__attr('T_BASELINE_END', baseline[1]))

        self.reject_criteria = reject_criteria

        self.epochs = mne.Epochs(
            self.raw, self.events, event_id=self.event_dict,
            tmin=tmin, tmax=tmax,
            reject=self.reject_criteria,
            baseline=baseline,
            reject_by_annotation=reject_by_annotation,
            preload=preload,
            **kwargs
        )

        # Visualise the log of dropped epochs
        if plot_drop_log:
            self.epochs.plot_drop_log()

        return self.epochs

    def block_average(self, rename=True):
        """Block averaging across trials."""
        # Dictionary with '<num_targets>/<hbo|hbr>' as keys and mne.Evoked object as value
        self.evoked_dict = {f'{event}/{ch_type}': self.epochs[event].average(picks=ch_type)
            for event in self.event_dict
            for ch_type in constants.HB_CHANNEL_TYPES
        }

        if rename:
            # Rename channels until the encoding of frequency in ch_name is fixed
            for condition in self.evoked_dict:
                self.evoked_dict[condition].rename_channels(lambda x: x.split(' ')[0])

        return self.evoked_dict

    def copy(self):
        """Get copy of member mne.raw."""
        return self.raw.copy()

    # @wrap
    def remove_backlight(raw, raw_backlight):
        """Backlight removal based on interpolation/smoothing."""
        raw = raw.copy()

        # Create design matrix of times (3rd order)
        regressors = sm.tools.tools.add_constant(np.c_[(times := raw.times), times**2, times**3]) # Timestamp (^1, ^2, ^3)

        # Fit RLM for every channel (row)
        fitted_backlight = np.apply_along_axis(
            lambda raw_backlight_ch: sm.RLM(raw_backlight_ch, regressors).fit().fittedvalues,
            1, raw_backlight
        )

        # Subtract predicted backlight signal from raw data of all wavelengths to remove backlight
        raw._data = raw.get_data() - np.repeat(fitted_backlight, int(len(raw.ch_names)/len(fitted_backlight)), axis=0)

        return raw

    def filter(self, l_freq=constants.F_L, h_freq=constants.F_H, l_trans_bandwidth=constants.L_TRANS_BANDWIDTH, h_trans_bandwidth=constants.H_TRANS_BANDWIDTH, filter_length='auto', **kwargs):
        return mne.filter.FilterMixin.filter(
            self.raw,
            l_freq=self.__attr('l_freq', l_freq),
            h_freq=self.__attr('h_freq', h_freq),
            l_trans_bandwidth=self.__attr('l_trans_bandwidth', l_trans_bandwidth),
            h_trans_bandwidth=self.__attr('h_trans_bandwidth', h_trans_bandwidth),
            filter_length=filter_length if filter_length else self.shape[-1]
        )

    def save_short_channels(self, max_dist=constants.SS_MAX_DIST):
        """Initialize member storing short channel mne.raw instance."""
        if not max_dist:
            self.raw_ss = self.raw.copy().drop_channels(utils.find_long_channels(self.raw.ch_names)[0])
        else:
            self.raw_ss = mne_nirs.channels.get_short_channels(self.raw, max_dist=max_dist)

    def scalp_coupling_index(raw, threshold=constants.THRESHOLD_SCI, *, plot_sci_drops=False):
        """Pick only channels with scalp coupling index above given threshold."""
        raw = raw.copy()

        sci = mne.preprocessing.nirs.scalp_coupling_index(raw)
        raw.info['bads'] += [ch_name for ch_name in raw.ch_names if ch_name in itertools.compress(raw.ch_names, sci < threshold)]

        if plot_sci_drops:
            plt.hist(sci)
            plt.xlabel = 'Scalp Coupling Index'
            plt.ylabel = 'Count'
            plt.xlim = [0, 1]

        return raw

    # @staticmethod
    def autopick_channels(raw, l_heart_rate=constants.L_HEART_RATE, h_heart_rate=constants.H_HEART_RATE, threshold_heart_rate=constants.THRESHOLD_HEART_RATE,
                          n_fft=None, ma_size=constants.MA_SIZE, *,
                          preserve_pairs=True, show_discarded=False, show_failed=False):
        """Automatic channel selection -- Heart-rate based.
        # > Fit Gaussian curve on the frequency spectrum of *HbO* between 0.6 and 1.8 Hz and filter out signals with low signal power (0.12 dB).
        # > Perdue, K. L.,Westerlund, A.,McCormick, S. A., and Nelson, C. A. (2014).
        # > Extraction of heart rate from functional near-infrared spectroscopy in infants.
        # > Journal of Biomedical Optics 19, 067010. doi:10.1117/1.JBO.19.6.067010
        # > https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4073682/
        # 
        # TODO: Automatic channel selection -- RMS-threshold based.

        Returns
        -------
        list
            Picks.
        """
        if isinstance(raw, NIRS):
            raw = raw.raw

        psd = raw.compute_psd(n_fft=n_fft if n_fft else len(raw))

        cut = np.argwhere((psd.freqs >= l_heart_rate) & (psd.freqs <= h_heart_rate)).squeeze()
        f = psd.freqs[cut]

        # Gaussian + constant
        _offset_gaussian = lambda x, a, x0, sigma, b: a * np.exp(-(x - x0)**2 / 2 / sigma**2) + b

        discards = set()
        failed = set()
        p0s = dict()
        popts = dict()
        for ch, ch_name in enumerate(psd.ch_names):
            if utils.is_long_channel(ch_name):
                # Calculate the log of the data (base 10)
                y = np.log(np.e) * np.log(psd.get_data()[ch][cut])

                # Smooth data using moving average
                y = sc.ndimage.uniform_filter1d(y, size=ma_size if ma_size else round(np.log(len(y))))

                # Intial values
                f0 = f[np.argmax(y)]
                b0 = np.median(y)
                a0 = np.max(y) - b0
                sigma_y_0 = np.median(np.abs(y - np.mean(y)))
                sigma_f_0 = (y - b0)*(y > b0) @ np.abs(f - f0) / np.sum(np.abs(y - b0))

                p0s[ch] = (a0, f0, sigma_f_0, b0)

                # Bound of parameters
                lower_bound = (0, f[0], 0, b0 - sigma_y_0)
                upper_bound = (np.inf, f[-1], np.ptp(f)/2, b0 + sigma_y_0)

                try:
                    popts[ch], pcov = sc.optimize.curve_fit(_offset_gaussian, f, y, p0s[ch], bounds=(lower_bound, upper_bound))
                except RuntimeError:
                    logging.warn(f'''Could not fit Gaussian curve for channel {psd.ch_names[ch]}. Discarding...''')
                    discards.add(ch)
                    failed.add(ch)
                else:
                    if popts[ch][0] < threshold_heart_rate:
                        discards.add(ch)

        if discards:
            match preserve_pairs:
                case False:
                    discards = discards.union(utils.find_ch_pairs(psd.ch_names, discards))
                case True:
                    discards -= set(utils.find_ch_unpaired(discards, psd.ch_names))
                case None:
                    pass

        def _make_overlay_plots(channels, title):
            fig, axs = plt.subplots(np.ceil(len(channels)/3).astype(int), min(3, len(channels)), figsize=(6 * min(3, len(channels)), 3 * np.ceil(len(channels)/3)), sharex=True, sharey=True)
            for ax, ch in zip(np.atleast_2d(axs).ravel(), channels):
                ax.plot(f, np.log(np.e) * np.log(psd.get_data()[ch][cut]), '+:b', markersize=3, alpha=0.6)
                ax.plot(f, _offset_gaussian(f, *p0s[ch]), 'o:g', markersize=2, label="$P_0$")
                if ch in popts:
                    ax.plot(f, _offset_gaussian(f, *popts[ch]), 'x:r', markersize=2, label="$P_{opt}$")
                    ax.set_title(f"{psd.ch_names[ch]} | $a = {popts[ch][0]:.2f}$")
                else:
                    ax.set_title(f"{psd.ch_names[ch]} | $a = N/A$")
                ax.legend()
            # fig.subplots_adjust(top=0.88)
            plt.suptitle(f"{title} | {len(channels)}")
            plt.tight_layout()

        if show_discarded and discards:
            _make_overlay_plots(discards, "Discarded Fits")

        if show_failed and failed:
            _make_overlay_plots(failed, "Failed Fits")

        raw.info['bads'] += [ch_name for ch, ch_name in enumerate(psd.ch_names) if ch in discards]

        if len(discards) == len(utils.find_long_channels(psd.ch_names)[0]):
            raise ValueError("No channels with heart beat signal power above threshold. Either decrease threshold or investigate using plots.")
        else:
            return [psd.ch_names[ch] for ch in set(range(len(psd.ch_names))) - discards]

    def default_pipeline(
            self,
            savepoints=dict(),
            remove_backlight=True,
            tddr=True,
            scalp_coupling_index=True,
            autopick_channels=True,
            short_channel_regression=True,
            pick_long_channels=True,
            bandpass=True,
            negative_correlation_enhancement=True,
            threshold_sci=constants.THRESHOLD_SCI,
            l_heart_rate=constants.L_HEART_RATE,
            h_heart_rate=constants.H_HEART_RATE,
            threshold_heart_rate=constants.THRESHOLD_HEART_RATE,
            n_fft=None,
            ma_size=constants.MA_SIZE,
            preserve_pairs=True,
            show_discarded=False,
            show_failed=False,
            ppf=constants.PPF,
            l_freq=constants.F_L,
            h_freq=constants.F_H,
            l_trans_bandwidth=constants.L_TRANS_BANDWIDTH,
            h_trans_bandwidth=constants.H_TRANS_BANDWIDTH,
            filter_length='auto',
            **kwargs
        ):
        """Default pipeline that runs a bunch of typical pre-processing functions and returns intermediate mne.raw instances as a dictionary.
        Stages: CW     (raw signal)
                CWx    (backlight removed raw signal)
                OD     (optical density)
                TDDR   (motion artifact removal)
                SCI    (scalp coupling index)
                AP     (autopick channels)
                SSR    (short-channel regression)
                HB     (chromophore/haemoglobin)
                LS     (pick long channels, after saving short channels)
                FL     (bandpass filtering)
                NCE    (negative correlation improvement)
        """
        if any(ch_type != 'fnirs_cw_amplitude' for ch_type in self.raw.info.get_channel_types()):
            raise ValueError("The default pipeline works only with channels of type fnirs_cw_amplitude.")
        self.process(
            # Save raw (CW amplitude) signals
                NIRS.save(savepoints)('CW'),
            # Remove Backlight
                NIRS.wrap(NIRS.remove_backlight)(self.raw_backlight, execute=remove_backlight),
                NIRS.save(savepoints)('CWx'),
            # Convert raw (CW amplitude) to optical density (OD) signals
                NIRS.wrap(mne.preprocessing.nirs.optical_density)(),
                NIRS.save(savepoints)('OD'),
            # Motion artifact removal -- Temporal Derivative Distribution Repair (TDDR)
                NIRS.wrap(mne.preprocessing.nirs.tddr)(execute=tddr),
                NIRS.save(savepoints)('TDDR'),
            # Pick only channels with high enough scalp coupling index
                NIRS.wrap(NIRS.scalp_coupling_index)(threshold_sci, execute=scalp_coupling_index),
                NIRS.save(savepoints)('SCI'),
            # Pick only channels with enough heart rate signal
                NIRS.wrap(NIRS.autopick_channels)(l_heart_rate, h_heart_rate, threshold_heart_rate, n_fft, ma_size, preserve_pairs=preserve_pairs,
                                                  show_discarded=show_discarded, show_failed=show_failed, execute=autopick_channels),
                NIRS.save(savepoints)('AP'),
            # Short-channel regression
                NIRS.wrap(mne_nirs.signal_enhancement.short_channel_regression)(max_dist=constants.SS_MAX_DIST, execute=short_channel_regression),
                NIRS.save(savepoints)('SSR'),
            # Optical Densities -> HbO and HbR concentrations -- Modified Beer Lambert Law (MBLL)
                # NIRS.wrap(mne.preprocessing.nirs.beer_lambert_law, ppf=0.1),
                NIRS.wrap(mbll.modified_beer_lambert_law)(ppf=self.__attr('PPF', ppf)),
                NIRS.save(savepoints)('HB'),
            # Pick long channels
                # Picking long channels removes all short channels, so before moving to that step, the short channels must be preserved
                NIRS.save_short_channels,
                NIRS.wrap(mne_nirs.channels.get_long_channels)(min_dist=constants.SS_MAX_DIST, max_dist=constants.LS_MAX_DIST, execute=pick_long_channels),
                NIRS.save(savepoints)('LS'),
            # Filter frequencies outside hemodynamic response range
                NIRS.wrap(mne.filter.FilterMixin.filter)(
                    l_freq=l_freq,
                    h_freq=h_freq,
                    l_trans_bandwidth=l_trans_bandwidth,
                    h_trans_bandwidth=h_trans_bandwidth,
                    filter_length=filter_length if filter_length else self.shape[-1],
                    execute=bandpass
                ),
                NIRS.save(savepoints)('FL'),
            # Negative correlation enhancement
                NIRS.wrap(mne_nirs.signal_enhancement.enhance_negative_correlation)(execute=negative_correlation_enhancement),
                NIRS.save(savepoints)('NCE')
        )
        return savepoints

    def get_psd(self, n_fft=None, ma_size=constants.MA_SIZE, **kwargs):
        if n_fft is None:
            n_fft = len(self.raw)

        psd = self.raw.compute_psd(n_fft=n_fft, **kwargs)

        if ma_size is None:
            ma_size = round(np.log(n_fft))

        psd._data = sc.ndimage.uniform_filter1d(psd.get_data(), size=ma_size)

        return psd

    psd = property(get_psd)

    def plot(self, duration=None, **kwargs):
        """Plot raw signals."""
        if duration is None:
            duration = self.DUR['exp']/3

        self.raw.plot(show_scrollbars=False, duration=duration, **kwargs)

    def plot_psd(self, n_fft=None, ma_size=constants.MA_SIZE, average=False, title="", **kwargs):
        """View power spectral densities of the signals."""
        fig = self.get_psd(n_fft, ma_size).plot(average=average, **kwargs)
        fig.suptitle(title)
        fig.subplots_adjust(top=0.88)

        return fig

    def plot_events(self, **kwargs):
        """Plot events in a scatter plot."""
        fig = mne.viz.plot_events(self.events, event_id=self.event_dict, sfreq=self.raw.info['sfreq'], show=False)
        fig.axes[0].set_yticklabels(self.cases)
        fig.subplots_adjust(right=0.7)

        return fig

    def plot_boxcar(self, title='', *, fig=None, axs=None, **kwargs):
        """Plot events in a boxcar plot."""
        if (fig is None) or (axs is None):
            fig, ax = plt.subplots(1, 1, figsize=(15, 6))

        plt.plot(self.raw.times, mne_nirs.experimental_design.create_boxcar(self.raw), axes=ax)
        plt.xlabel("Time (s)")
        plt.title(title)
        plt.legend(self.cases, loc='upper right')

        return fig

    def plot_sensors_3d(self, **kwargs):
        """Show sensors on fsaverage brain."""
        subjects_dir = os.path.join(mne.datasets.sample.data_path(), 'subjects')
        mne.datasets.fetch_fsaverage(subjects_dir=subjects_dir)
        brain = mne.viz.Brain('fsaverage', subjects_dir=subjects_dir, alpha=0.5, cortex='low_contrast')
        brain.add_head()
        brain.add_sensors(self.raw.info, trans='fsaverage')
        brain.show_view(azimuth=90, elevation=90, distance=500)

    def plot_average_heatmap(self, picks=None, exclude='bads', clim={'hbo': [-10, 10], 'hbr': [-10, 10]}, *, fig=None, axs=None, **kwargs):
        """Plot heatmap of block averaged signals for all channels, for all cases."""
        if (fig is None) or (axs is None):
            fig, axs = plt.subplots(2, len(self.cases), figsize=(18, 6))

        for ax, event in zip(np.atleast_2d(axs.T), self.cases):
            self.epochs[event].average().plot_image(picks=picks, exclude=exclude, axes=ax, clim=clim, show=False)
            ax[0].set_xlabel(None)
            ax[0].set_title(f"{event} Targets | HbO")
            ax[1].set_title(f"{event} Targets | HbR")
            for ax_i in ax:
                ax_i.axvline(0, c='k', ls='--', lw=0.9)

        fig.suptitle("Block-Averaged Signals Across Trials for Channels and Number of Targets")
        return fig

    def plot_average_waveform(self, picks=None, cases=None, separate_cases=False, *, title_format_hex=True, fig=None, axs=None, sharex=True, sharey=True, **kwargs):
        """Plot block averaged signals for given picks and cases."""
        if cases is None:
            cases = self.cases

        s_d = utils.get_s_d(picks)

        if (fig is None) or (axs is None):
            if separate_cases:
                fig, axs = plt.subplots(len(s_d), len(cases), figsize=(6 * len(cases), 3 * len(s_d)), sharey=sharey, sharex=sharex)
            else:
                fig, axs = plt.subplots(np.ceil(len(s_d)/3).astype(int), min(3, len(s_d)),
                                        figsize=(6 * min(3, len(s_d)), 3 * np.ceil(len(s_d)/3)), sharey=sharey, sharex=sharex)   

        if separate_cases:
            for ax, case in zip(np.atleast_2d(axs.T), cases):
                evoked = self.epochs[case].average(picks=picks)
                for ax_i, s_d_i in zip(ax, s_d):
                    ax_i.axvline(0, color='k', linestyle='--', alpha=0.5)
                    for ch_type in constants.HB_CHANNEL_TYPES:
                        color = 'r' if ch_type == 'hbo' else 'b' if ch_type == 'hbr' else 'g'
                        if (ch_name := f'{s_d_i} {ch_type}') in picks:
                            ax_i.plot(evoked.times, evoked.get_data(picks=ch_name).squeeze().T * 1e6, color=color)
                            ax_i.plot(self.epochs.times, self.epochs[case].get_data(picks=ch_name).squeeze().T * 1e6, color=color, alpha=0.1)
                            ax_i.set_title(f"{case} Targets | {utils.dec_to_hex([s_d_i])[0] if title_format_hex else s_d_i}")
                            ax_i.set_xlabel("Times")
                            ax_i.set_ylabel(r"$\mu M$")
        else:
            for ax in axs.ravel():
                ax.axvline(0, color='k', linestyle='--', alpha=0.5)
            for i, case in enumerate(cases):
                evoked = self.epochs[case].average(picks=picks)
                for ax_i, s_d_i in zip(axs.ravel(), s_d):
                    for ch_type in constants.HB_CHANNEL_TYPES:
                        color = 'r' if ch_type == 'hbo' else 'b' if ch_type == 'hbr' else 'g'
                        if (ch_name := f'{s_d_i} {ch_type}') in picks:
                            alpha = (i + 1)/(len(cases) + 1)
                            label = case if ch_type == 'hbo' else None
                            ax_i.plot(evoked.times, evoked.get_data(picks=ch_name).squeeze().T * 1e6, color=color, alpha=alpha, label=label)
                            # ax_i.plot(self.epochs.times, self.epochs[case].get_data(picks=ch_name).squeeze().T * 1e6, color=color, alpha=0.2)
                            ax_i.set_title(f"{utils.dec_to_hex([s_d_i])[0] if title_format_hex else s_d_i}")
                            ax_i.legend(loc='upper right')
                            ax_i.set_xlabel("Times")
                            ax_i.set_ylabel(r"$\mu M$")

        return fig

if __name__ == '__main__':
    pass
