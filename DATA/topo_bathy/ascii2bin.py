#!/usr/bin/env python

from __future__ import print_function
import sys

try:
    import numpy as np
except ImportError:
    print('This script requires NumPy.')
    sys.exit()

if '--help' in sys.argv or '-h' in sys.argv or len(sys.argv) != 3:
    print('Usage: %s <input> <output>' % (sys.argv[0],))
    sys.exit()

# Read input file
print('Reading input file %s ...' % (sys.argv[1],))
data = np.genfromtxt(sys.argv[1])

# Convert to 16-bit integers
data2 = data.astype(np.int16)
if any(data != data2):
    print('Warning: Data set does not fit in signed 16-bit integers!')

# Add a byte-order mark
byteorder = np.array([0x1234], dtype=np.int16)
data2 = np.concatenate((byteorder, data2))

# Save output file
print('Writing output to file %s ...' % (sys.argv[2],))
data2.tofile(sys.argv[2])

print('Done!')

