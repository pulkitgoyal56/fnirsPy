# Base directory of data
DATA_DIR = 'data'

# Project name
PROJECT = 'Multi-object Tracking'

class DEVICE(object):
    # Number of wavelengths from which data is collected
    N_WAVELENGTHS = 4

    # Maximum length of a short channel
    SS_MAX_DIST = 0.01 # m

    # Maximum length of a long channel
    LS_MAX_DIST = 0.10 # m

    # Time drift factor
    TIME_DRIFT_FACTOR = 1.0045 

    # Additional device info
    INFO = {
        'type': 'fNIRS-CW',
        'model': 'optoHIVE'
    }
    EXPERIMENTER = 'optoHIVE Team'


# Maximum threshold for channel rejection
# Reject epochs based on maximum peak-to-peak (PTP) signal amplitude, i.e. the absolute difference between the lowest and the highest signal value.
# If the PTP signal amplitude of any one channel exceeds the rejection threshold, the respective epoch will be dropped.
REJECT_CRITERIA = {
    'hbo': 80e-6
}

F_L = 0.02
F_H = 0.4
L_TRANS_BANDWIDTH = 0.01
H_TRANS_BANDWIDTH = 0.2

HB_CHANNEL_TYPES = ['hbo', 'hbr']

# The partial pathlength factors for different frequencies
PPF = {
    770: 6.18052629,
    774: 6.18052629,
    810: 5.930508,
    817: 5.930508,
    855: 5.50968514,
    865: 5.50968514,
    885: 5.17117029,
    892: 5.17117029
}


# Expected locations of reference points (at least 4 points are required to recreate mapping)
MNI_REFERENCE_LOCATIONS = {
    'nasion': (0.0083, 86.8110, -39.9830),
    'lpa': (-86.0761, -19.9897, -47.9860),
    'rpa': (85.7939, -20.0093, -48.0310),
    'Cz': (0.4009, -9.1670, 100.2440)
}

DEFAULT_REFERENCE_LOCATIONS = MNI_REFERENCE_LOCATIONS
DEFAULT_REFERENCE = 'mri'
