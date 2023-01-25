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

        self.N_WAVELENGTHS = self.DEVICE.N_WAVELENGTHS
        # self.LS_MAX_DIST = self.DEVICE.LS_MAX_DIST
        # self.SS_MAX_DIST = self.DEVICE.SS_MAX_DIST
        self.TIME_DRIFT_FACTOR = self.DEVICE.TIME_DRIFT_FACTOR

    def __attr(self, attr, var):
        if var is not None:
            self.__setattr__(attr, var)
        return self.__getattribute__(attr)

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

    def pick_wavelengths(self, wavelengths_picked=None):
        wavelengths_picked = self.__attr('WAVELENGTHS_PICKED', wavelengths_picked)

        # Indices of all the channels available (beware, these are not the same as the initial channel numbers!)
        # picks = mne.pick_types(self.raw.info, meg=False, fnirs=True) # Select channels with picked wavelengths

        # Pick long channels (for picked wavelength)
        # self.raw.pick([ch for ch in utils.find_long_channels(self.raw.ch_names)[0] if int(ch.split()[1]) in self.WAVELENGTHS_PICKED])

        self.raw.pick([ch for ch in self.raw.ch_names if int(ch.split()[1]) in wavelengths_picked])

    def set_bad(self, bad_channels):
        self.BAD_CHANNELS = bad_channels

        if isinstance(bad_channels[0], int):
            bad_channels = [self.CH_NAMES[ch] for ch in bad_channels]

        self.raw.info['bads'] = bad_channels

    def read_config(self, config_file_path):
        """Read additional configuration."""
        self.config_file_path = pathlib.Path(config_file_path).with_suffix('.toml')

        with open(self.config_file_path, 'rb') as f:
            self.config = tomli.load(f)

            self.WAVELENGTHS_PICKED = self.config['WAVELENGTHS_PICKED']
            self.PPF = [float(constants.PPF[wavelength]) for wavelength in self.WAVELENGTHS_PICKED]
            self.N_PROBES = int(self.config['N_PROBES'])
            self.N_HEMISPHERES = int(self.config['N_HEMISPHERES'])
            self.S_D = utils.hex_to_dec(self.config['S_D'])
            self.CH_UNUSED = set(self.config['CH_UNUSED'])
            self.T_EXP_START = float(self.config['T_EXP_START'])
            self.T_EPOCH_START = float(self.config['T_EPOCH_START'])
            self.T_EPOCH_END = float(self.config['T_EPOCH_END'])
            self.T_BASELINE_START = float(self.config['T_BASELINE_START'])
            self.T_BASELINE_END = float(self.config['T_BASELINE_END'])

    def read_raw_fif(self, raw_file_path, config_file_path=None, pick_wavelengths=True, remove_backlight=True):
        self.raw_file_path = raw_file_path.parent / pathlib.Path(raw_file_path.stem.split('.')[0]).with_suffix('.raw.fif')

        if config_file_path:
            self.read_config(config_file_path)

        # Read '.raw.fif' file
        raw_fif = mne.io.read_raw_fif(self.raw_file_path, preload=True)
        # raw_fif.crop(tmin=120) # Delete first 60s for this dataset (idle data)

        # Wavelengths available
        self.WAVELENGTHS = pd.unique([int(ch.split()[1]) for ch in raw_fif.ch_names])

        # Total number of wavelengths
        self.N_WAVELENGTHS_T = len(self.WAVELENGTHS)

        # Number of wavelengths used in analysis
        self.N_WAVELENGTHS = len(self.WAVELENGTHS_PICKED) # 2

        # # Maximum number of channels
        # self.M_CHANNELS = int(len(raw_fif.ch_names) / self.N_WAVELENGTHS_T)  # per wavelength # self.N_PROBES ** 2

        # # Number of probes
        # self.M_PROBES = np.ceil(np.sqrt(self.M_CHANNELS) / self.N_HEMISPHERES)  # per hemisphere # self.N_PROBES

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
        raw_fif.drop_channels([raw_fif.ch_names[ch] for ch in self.CH_UNUSED])

        # Names of the wavelength specific channels (only used channels)
        self.CH_NAMES = raw_fif.ch_names

        # Sampling frequency (based on difference between timestamps in consecutive readings ~54ms)
        self.F_S = raw_fif.info['sfreq'] # * self.TIME_DRIFT_FACTOR          # fNIRS recording frequency, in Hertz

        # Set recording start and end times
        self.T_REC_START = 0                                                 # fNIRS recording start time, in seconds
        self.T_REC_END = len(raw_fif)/self.F_S                               # fNIRS recording end time, in seconds

        # Read '-backlight.raw.fif' file
        # Backlight intensities (for used channels only)
        if remove_backlight:
            backlight_file_path = raw_file_path.parent / pathlib.Path(raw_file_path.stem.split('.')[0] + '-backlight').with_suffix('.raw.fif')
            self.raw_backlight = mne.io.read_raw_fif(backlight_file_path, preload=True).get_data()

        self.raw = raw_fif

        if pick_wavelengths:
            self.pick_wavelengths()

        return self.raw

    def read_raw_csv(self, raw_file_path, config_file_path, pick_wavelengths=True, remove_backlight=True):
        self.raw_file_path = pathlib.Path(raw_file_path).with_suffix('.csv')

        if config_file_path:
            self.read_config(config_file_path)

        # Read CSV data as Pandas DataFrame
        data_pd = pd.read_csv(self.raw_file_path)

        # Wavelengths available (automatic extraction)
        self.WAVELENGTHS = [int(match.groups()[0]) for column in data_pd.columns if (match := re.compile(r'(\d+)\[nm\]').match(column))] # [855, 770, 810, 885]

        # Total number of wavelengths
        self.N_WAVELENGTHS_T = len(self.WAVELENGTHS)

        # Number of wavelengths used in analysis
        self.N_WAVELENGTHS = len(self.WAVELENGTHS_PICKED) # 2

        # # Number of probes
        # self.M_PROBES = self.N_PROBES # per hemisphere

        # # Maximum number of channels
        # self.M_CHANNELS = (self.M_PROBES)**2 * self.N_HEMISPHERES # per wavelength

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

        # Remove unused channels and create a new DataFrame
        data_pd = data_pd.loc[data_pd['Channel'].isin(self.CH_USED)]

        # # Check if the number of channels are as expected (they must not be more than `N_CHANNELS`)
        # if len(data_pd['Channel'].unique()) > M_CHANNELS:
        #     raise ValueError(f'''Duplicate channels. Expected (max.) - {M_CHANNELS}; Received - {len(data_pd["Channel"].unique())}.\n
        #                          Please pick one of the duplicates and mark the others in `CH_UNUSED` in the related config file.''')

        # Names of the wavelength specific channels
        # N_CHANNELS * N_WAVELENGTHS_T
        self.CH_NAMES = [f'{self.CH_MAP[ch]} {wavelength}'
                    for ch in list(self.CH_MAP.keys()) if ch not in self.CH_UNUSED
                    for wavelength in self.WAVELENGTHS]

        # Set recording start and end times
        self.T_REC_START = -data_pd['Time[ms]'].iloc[0]/1000/self.TIME_DRIFT_FACTOR    # fNIRS recording start time, in seconds
        self.T_REC_END = np.ptp(data_pd['Time[ms]'])/1000/self.TIME_DRIFT_FACTOR       # fNIRS recording end time, in seconds

        # Sampling frequency (based on difference between timestamps in consecutive readings ~54ms)
        self.F_S = (len(data_pd) - len(self.S_D))/len(self.S_D)/self.T_REC_END         # fNIRS recording frequency, in Hertz

        ## Create mne.Info object
        # Type of channels -- Raw fNIRS Continuous Wave Amplitude
        # https://mne.tools/stable/glossary.html#term-data-channels
        CH_TYPES = 'fnirs_cw_amplitude' # [] * N_CHANNELS * N_WAVELENGTHS_T * N_HEMISPHERES

        # Create MNE.Info Object
        config_csv = mne.create_info(ch_names=self.CH_NAMES, sfreq=self.F_S, ch_types=CH_TYPES)

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

        ## Create mne.io.Raw object
        raw_csv = mne.io.RawArray(data_np, config_csv)

        # Backlight intensities (for used channels only)
        if remove_backlight:
            self.raw_backlight = data_pd['BL'].to_numpy().reshape(-1, self.N_CHANNELS).T

        self.raw = raw_csv

        if pick_wavelengths:
            self.pick_wavelengths()

        return self.raw

    def read_raw(self, raw_file_path, config_file_path=None, pick_wavelengths=True, remove_backlight=True):
        raw_file_path = pathlib.Path(raw_file_path)

        if raw_file_path.suffix == '':
            # Check if file with file name exists, with any extention
            if (raw_file_path := list(raw_file_path.parent.glob(f'{raw_file_path.stem}.*'))):
                raw_file_path = raw_file_path[0]
            else:
                raise FileNotFoundError(f'No file of the sort `{raw_file_path}.*`')

        match raw_file_path.suffix:
            case '.csv':
                return self.read_raw_csv(raw_file_path, config_file_path, pick_wavelengths, remove_backlight)
            case '.fif':
                return self.read_raw_fif(raw_file_path, config_file_path, pick_wavelengths, remove_backlight)
            case other:
                raise ValueError(f'Unsupported fNIRS file format - {other}')

    def read_annotation(self, annotation_file_path):
        self.annotation_file_path = pathlib.Path(annotation_file_path).with_suffix('.mat')

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
        self.DUR = self.mat.filter(regex=(".*_p$")).mean().round().rename(lambda c_n: c_n[:-2]) # .astype(int)
        self.DUR['trial'] = sum(self.DUR)

        # Read experiment end time, relative to its start time, vis-à-vis its duration
        endtime_file_path = annotation_file_path.parent / pathlib.Path(annotation_file_path.stem.rsplit('_', 1)[0] + '_endtime').with_suffix('.mat')
        self.DUR['exp'] = float(sc.io.loadmat(endtime_file_path)['expEnd'])  # <expEnd - exp> -- duration of the entire experiment, in seconds

        # Set the duration of the recording
        self.DUR['rec'] = self.T_REC_END                                     # recording duration, in seconds

        # Read experiment end time and set start time
        # self.T_EXP_START = 0               # <exp>                         # experiment start time, in seconds; offset due to trigger delay, in seconds
        self.T_EXP_END = self.T_EXP_START + self.DUR['exp']                  # experiment end time, in seconds

        # Set annotations in the Raw object
        self.raw.set_annotations(mne.Annotations(
            onset=self.T_EXP_START + self.mat['motion_e'] - self.T_REC_START,
            duration=[self.DUR['motion']] * len(self.mat),
            description=self.mat['num_targets'].astype(int)
        ))

    def read_montage(self, montage_file_path):
        self.montage_file_path = pathlib.Path(montage_file_path).with_suffix('.elc')
        self.raw.set_montage(mne.channels.read_custom_montage(self.montage_file_path, coord_frame='mri'))

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

        return self

    def process(self, *funcs):
        if len(funcs) > 1:
            for func in funcs:
                self.process(func)
        else:
            if isinstance((raw := funcs[0](self)), type(self.raw)):
                self.raw = raw

    @staticmethod
    def save(savepoints):
        """Saves objects passed in a dictionary."""
        def wrapper(label):
            def subwrapper(self):
                savepoints[label] = self.raw
            return subwrapper
        return wrapper

    @staticmethod
    def wrap(func, *args, **kwargs):
        def wrapper(self):
            return func(self.raw, *args, **kwargs)
        return wrapper

    def get_epochs(self, tmin=None, tmax=None, baseline=(None, None), reject_criteria=constants.REJECT_CRITERIA, reject_by_annotation=False, **kwargs):
        """Extract Epochs."""

        tmin = self.__attr('T_EPOCH_START', tmin)
        tmax = self.__attr('T_EPOCH_END', tmax)

        baseline = (self.__attr('T_BASELINE_START', baseline[0]), self.__attr('T_BASELINE_END', baseline[1]))

        self.reject_criteria = reject_criteria

        # Extract events of interest
        self.events, self.event_dict = mne.events_from_annotations(self.raw)

        self.epochs = mne.Epochs(
            self.raw, self.events, event_id=self.event_dict,
            tmin=self.T_EPOCH_START, tmax=self.T_EPOCH_END,
            reject=reject_criteria,
            baseline=baseline,
            reject_by_annotation=reject_by_annotation,
            preload=True,
            **kwargs
        )

        ## Visualise the log of dropped epochs
        # epochs.plot_drop_log()

        return self.events, self.event_dict, self.epochs

    def block_average(self, rename=True):
        # Block averaging across trials
        # Dictionary with '<num_targets>/<hbo|hbr>' as keys and mne.Evoked object as value
        self.evoked_dict = {f'{event}/{ch_type}': self.epochs[event].average(picks=ch_type)
            for event in self.event_dict.keys()
            for ch_type in constants.HB_CHANNEL_TYPES
        }

        if rename:
            # Rename channels until the encoding of frequency in ch_name is fixed
            for condition in self.evoked_dict:
                self.evoked_dict[condition].rename_channels(lambda x: x[:-4])

        return self.evoked_dict

if __name__ == '__main__':
    pass
