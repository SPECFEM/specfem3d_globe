#!/usr/bin/env python

# conversion of topography files from ASCII to binary to reduce their size.
# endianness (little-endian or big-endian machines) is properly handled based on a byte-order mark written in the header.
# script written by Elliott Sales de Andrade.

from __future__ import print_function
import sys
from io import BytesIO

try:
    import numpy as np
except ImportError:
    print('This script requires NumPy.')
    sys.exit()

if '--help' in sys.argv or '-h' in sys.argv or len(sys.argv) != 3:
    print('Usage: %s <input> <output>' % (sys.argv[0],))
    sys.exit()

# Output config
print('Reading input file %s ...' % (sys.argv[1],))
print('Writing output to file %s ...' % (sys.argv[2],))

# Input file
inf = open(sys.argv[1], 'rb')

# Output file
outf = open(sys.argv[2], 'wb')

# Add a byte-order mark
byteorder = np.array([0x1234], dtype=np.int16)
byteorder.tofile(outf)

i = 0
while True:
    # Read input file
    text = inf.readlines(1024*1024)  # 1M at a time
    if not len(text):
        break
    ioin = BytesIO(''.join(text))
    data = np.genfromtxt(ioin)
    if not len(data):
        break

    # Convert to 16-bit integers
    data2 = data.astype(np.int16)
    if any(data != data2):
        print('Warning: Data set does not fit in signed 16-bit integers!')

    # Save output file
    data2.tofile(outf)

    i = i + 1
    if i % 10 == 0:
        print('%d MB ...' % (i,))

inf.close()
outf.close()
print('Done!')

