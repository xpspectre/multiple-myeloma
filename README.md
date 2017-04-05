# Multiple Myeloma Project

Note: I used Python 3 for this.

## Setup

Make a new Anaconda3 environment
```
conda create -n mm-env python=3.5 numpy scipy matplotlib pandas msgpack-python scikit-learn pillow xlrd
source activate mm-env
pip install tensorflow msgpack-numpy
```

Make a symlink in this directory to the data store
```
ln -s data <mm data dir>
```