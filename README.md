# Multiple Myeloma Project

## Setup

Make a new Anaconda3 environment
```
conda create -n mm-env numpy scipy matplotlib pandas msgpack-python scikit-learn pillow
source activate mm-env
pip install tensorflow
```

Make a symlink in this directory to the data store
```
ln -s <mm data dir> data
```