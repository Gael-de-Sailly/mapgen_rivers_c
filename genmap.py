#!/usr/bin/env python3

import numpy as np
import noise
import sys

if len(sys.argv) > 1:
    mapsize = int(sys.argv[1])
else:
    mapsize = 400

scale = mapsize / 2
n = np.zeros((mapsize+1, mapsize+1))

# Set noise parameters
params = {
    "octaves" : int(np.ceil(np.log2(mapsize)))+1,
    "persistence" : 0.5,
    "lacunarity" : 2.,
}

# Determine noise offset randomly
xbase = np.random.randint(1024)
ybase = np.random.randint(1024)

# Generate the noise
for x in range(mapsize+1):
    for y in range(mapsize+1):
        n[x,y] = noise.snoise2(x/scale + xbase, y/scale + ybase, **params)

nn = n*mapsize/5 + mapsize/20

filename = "dem"
if len(sys.argv) > 2:
    filename = sys.argv[2]

with open(filename, 'wb') as f:
    f.write(np.array(nn.shape, dtype='<u2').tobytes())
    f.write(nn.astype('<f8').tobytes())


print(nn)
