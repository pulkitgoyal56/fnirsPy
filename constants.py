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
    TIME_DRIFT_FACTOR = 1.01 

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

HB_CHANNEL_TYPES = ['hbo', 'hbr']
