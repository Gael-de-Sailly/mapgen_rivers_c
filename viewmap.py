#!/usr/bin/env python3

import numpy as np
import noise
import sys
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = 'dem'

default_types = {
    'dem': '<f8',
    'dem_new': '<f8',
    'dirs': '<u1',
    'rivers': '<f8',
    'lakes': '<f8',
}

if len(sys.argv) > 2 and sys.argv[2] != 'log':
    data_dtype = '<' + sys.argv[2]
else:
    if filename in default_types:
        data_dtype = default_types[filename]
    else:
        data_dtype = '<f8'

islog = len(sys.argv) > 2 and sys.argv[-1] == 'log'

with open(filename, 'rb') as f:
    size_raw = f.read(4)
    size = tuple(np.frombuffer(size_raw, dtype='<u2'))
    data = np.frombuffer(f.read(), dtype=data_dtype).reshape(size)
    if islog:
        data = np.log(data)
    plt.imshow(data)
    plt.colorbar()
    plt.show()
