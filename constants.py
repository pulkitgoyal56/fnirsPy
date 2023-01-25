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

    # Additional device info
    INFO = {
        'type': 'fNIRS-CW',
        'model': 'optoHIVE'
    }
    EXPERIMENTER = 'optoHIVE Team'
