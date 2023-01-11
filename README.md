# fNIRS-attentional-load
fNIRS Data Analysis for Multi-object Tracking Task

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
