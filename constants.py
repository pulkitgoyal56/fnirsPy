# Base directory of data
DATA_DIR = 'data'

# Project name
PROJECT = 'Multi-object Tracking'

class DEVICE(object):
    # Number of wavelengths from which data is collected
    N_WAVELENGTHS = 4

    # Short separation
    SS_SEPARATION = 0.007

    # Long separation
    LS_SEPARATION = 0.03

    # Time drift factor
    TIME_DRIFT_FACTOR = 1.0045

    # Additional device info
    INFO = {
        'type': 'fNIRS-CW',
        'model': 'optoHIVE'
    }
    EXPERIMENTER = 'optoHIVE Team'

# Maximum length of a short channel
SS_MAX_DIST = 0.01 # m

# Maximum length of a long channel
LS_MAX_DIST = 0.10 # m

# Maximum threshold for channel rejection
# Reject epochs based on maximum peak-to-peak (PTP) signal amplitude, i.e. the absolute difference between the lowest and the highest signal value.
# If the PTP signal amplitude of any one channel exceeds the rejection threshold, the respective epoch will be dropped.
REJECT_CRITERIA = None # {'hbo': 80e-6}

# See recommendations by Pinti et. al. 2014
# https://www.frontiersin.org/articles/10.3389/fnhum.2018.00505/full
F_L = 0.01
F_H = 0.09
# The transition band is chosen such that the order of the resulting filter is >1000.
L_TRANS_BANDWIDTH = 0.006
H_TRANS_BANDWIDTH = 0.006

HB_CHANNEL_TYPES = ['hbo', 'hbr']

# The partial pathlength factors for different frequencies
# See - https://github.com/mne-tools/mne-python/pull/9843
PVC = 60 # 1 # Partial Volume Correction
PPF = {
    770: 6.18052629 / PVC,
    774: 6.18052629 / PVC,
    810: 5.930508 / PVC,
    817: 5.930508 / PVC,
    855: 5.50968514 / PVC,
    865: 5.50968514 / PVC,
    885: 5.17117029 / PVC,
    892: 5.17117029 / PVC
}

# Expected locations of reference points (at least 4 points are required to recreate mapping)
MNI_REFERENCE_LOCATIONS = {
    'nasion': (0.0083, 86.8110, -39.9830),
    'lpa': (-86.0761, -19.9897, -47.9860),
    'rpa': (85.7939, -20.0093, -48.0310),
    'Cz': (0.4009, -9.1670, 100.2440),
    'IPS': (34., -54., 40.)
}

DEFAULT_REFERENCE_LOCATIONS = MNI_REFERENCE_LOCATIONS
DEFAULT_REFERENCE = 'mri'
