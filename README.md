# `fnirsPy`

_fnirsPy_ is an open-source wrapper built on [MNE-Python](https://github.com/mne-tools/mne-python) (Gramfort et al., 2014) and [MNE-NIRS](https://github.com/mne-tools/mne-nirs) (Luke et al., 2021) that makes it easy to deploy complete yet customizable end-to-end data processing pipelines for custom near-infrared spectroscopy (NIRS) data with only a few lines of code. It is an imperative-style high-level API designed to be Pythonic, which focuses on minimum redundancy in code, so it augments and does not replace MNE, the functions of which users can still utilize harmoniously. It also aims to implement a standardized procedure for handling fNIRS data to foster reproducibility in analysis and includes a default pipeline based on growing recommendations for standardized fNIRS data preprocessing steps in literature so that users can focus on the scientific question without much toil over writing code. The package also implements new methods that are not available in the underlying MNE software suite. Furthermore, a few additional plotting functions are included with fNIRS publication conventions in mind along with some handy utilities that provide frequently used functions.

> The library is currently in pre-alpha and few of its functions are specifically developed for the in-development OptoHive fNIRS device by RELab (ETH Zürich) (Wyser et al., 2017).  
> It is planned to be developed further to be able to read any task data following the [BIDS](https://bids.neuroimaging.io/) standard for neuroimaging data.

## Setup (for Debian-Based Linux)

### Update System Packages

``` sh
sudo apt update
sudo apt -y upgrade
```

### Install Node

``` sh
sudo apt install -y nodejs
```

### Install Graphics Dependencies

> These are probably more than necessary

``` sh
sudo apt install -y ffmpeg libsm6 libxext6 python-qt4 libgl1-mesa-glx xvfb
```

### Install Miniconda

``` sh
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_22.11.1-1-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### Create Conda Environment

``` sh
conda env create -f environment.yml
```

### Run Jupyter Lab

``` sh
jupyter lab --ip 0.0.0.0 --no-browser
```
