#!/usr/bin/env python

# conversion of topography files from ASCII to binary to reduce their size.
# endianness (little-endian or big-endian machines) is properly handled based on a byte-order mark written in the header.
# script written by Elliott Sales de Andrade.

from __future__ import print_function
import sys
import os
from io import BytesIO

try:
    import numpy as np
except ImportError:
    print('This script requires NumPy.')
    sys.exit()


def convert_etopo_ascii2bin(filename_in,filename_out):
    # Output config
    print('Reading input file %s ...' % (filename_in,))
    print('Writing output to file %s ...' % (filename_out,))

    # Input file - ascii format
    inf = open(filename_in, 'rb')

    size_b = os.path.getsize(filename_in)
    size_mb = size_b / 1024.0 / 1024.0
    print('input file size: %6.1f MB' % (size_mb))

    # Output file - binary format
    outf = open(filename_out, 'wb')

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
            print('%6d MB ...' % (i))

    inf.close()
    outf.close()

    print('Done!')


if __name__ == '__main__':
    # gets arguments
    if '--help' in sys.argv or '-h' in sys.argv or len(sys.argv) != 3:
        print('Usage: %s <input> <output>' % (sys.argv[0],))
        sys.exit()

    # Input file
    filename_in = sys.argv[1]

    # Output file
    filename_out = sys.argv[2]

    convert_etopo_ascii2bin(filename_in,filename_out)

