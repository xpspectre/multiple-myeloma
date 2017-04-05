# Compressed msgpack with numpy library
# Importing this module also automatically applies the needed msgpack_numpy.patch()

import msgpack
import msgpack_numpy
import gzip
import numpy as np

msgpack_numpy.patch()


def packb(data):
    return gzip.compress(msgpack.packb(data))


def unpackb(stream):
    return msgpack.unpackb(gzip.decompress(stream))


# Test code
if __name__ == '__main__':

    data = {
        'a': [1,2,3],
        'b': np.zeros((50,), dtype=np.float32)
    }

    # Test uncompressed msgpack
    filename = 'testfile'
    with open(filename, 'wb') as f:
        f.write(msgpack.packb(data))
    with open(filename, 'rb') as f:
        data_out = msgpack.unpackb(f.read())

    # Test compressed msgpack
    filename_compressed = 'testfile.gz'
    with open(filename_compressed, 'wb') as f:
        f.write(packb(data))
    with open(filename_compressed, 'rb') as f:
        data_out = unpackb(f.read())

    1
