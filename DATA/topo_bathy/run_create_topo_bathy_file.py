#!/usr/bin/env python
#
# script to download topography data and create a SPECFEM readable topo/bathy file
#
# required python modules:
#   - requests
#   - timeit
#   - zipfile
#   - numpy
#
# required external packages:
#   - GMT       - e.g. version 5.4.*
#
from __future__ import print_function

import os
import sys
import subprocess

import requests
import timeit
import zipfile

## numpy
try:
    import numpy as np
except ImportError:
    print('This script requires NumPy.')
    sys.exit()

## GMT
## http://gmt.soest.hawaii.edu
## install by: > sudo apt install gmt
## setup gmt functions by:
## > source /usr/share/gmt/tools/gmt_functions.sh
try:
    print("checking GMT version:")
    cmd = 'gmt --version'
    print("> ",cmd)
    version = subprocess.check_output(cmd, shell=True)
except:
    print("Error using `gmt`")
    print("install by: > sudo apt install gmt")
    sys.exit(1)
# avoid bytes string issues with strings like b'Hello', converts to text string
if isinstance(version, (bytes, bytearray)): version = version.decode("utf-8")
version = version.strip()
print("GMT version: %s" % (version))
# get version numbers for later (grdconvert command format changes between version 5.3 and 5.4)
elem = version.split(".")
gmt_major = int(elem[0])
gmt_minor = int(elem[1])
if gmt_major >= 5 and gmt_minor >= 3: print("version ok")
print("")


####################################################################
# User parameters

## ETOPO1
# see: https://www.ngdc.noaa.gov/mgg/global/global.html
#
# ETOPO1 is a 1 arc-minute global relief model of Earth's surface
# that integrates land topography and ocean bathymetry.
# Built from global and regional data sets, it is available in "Ice Surface" (top of Antarctic and Greenland ice sheets)
# and "Bedrock" (base of the ice sheets).
#
# The grid-registered is the authoritative registration.
# (The cell-registered is derived from the grid-registered, and the conversion produces slightly flattened relief.)
#
# with matlab scripts:
#   has slightly different number of grid points 21601 x 10801
#filename_web1 = 'ETOPO1_Ice_g_geotiff.tif'
#url_etopo1 = 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/georeferenced_tiff/ETOPO1_Ice_g_geotiff.zip'
#
# cell-registered binary version:
filename_web1 = 'etopo1_ice_c_i2.bin'
url_etopo1 = 'http://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/cell_registered/binary/etopo1_ice_c_i2.zip'

## ETOPO2
# see: https://www.ngdc.noaa.gov/mgg/global/etopo2.html
#
# ETOPO2v2 is available in both a downloadable cell-centered version,
# named ETOPO2v2c (pixel registered, where the cell boundaries are lines of even minutes of latitude and longitude,
# centered on intersections of lines of odd minutes of latitude and longitude)
# and a grid-centered, version, available via design-a-grid (with cell boundaries defined by lines of odd minutes of latitude
# and longitude, meaning that cells were centered on the integer multiples of 2 minutes [even minutes] of latitude and longitude).
#
# The cell-centered grid is the authoritative version.
# (The grid-centered version was derived from the cell-centered grid and the conversion produces slightly flattened relief.
#  This is a change from the original ETOPO2 (2001), which was available as grid-registered only.)
#
# netCDF format
#filename_web2 = 'ETOPO2v2c_f4.nc'
#url_etopo2 = 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2/ETOPO2v2-2006/ETOPO2v2c/netCDF/ETOPO2v2c_f4_netCDF.zip'
# Linux little-endian
filename_web2 = 'ETOPO2v2c_i2_LSB.bin'
url_etopo2 = 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2/ETOPO2v2-2006/ETOPO2v2c/raw_binary/ETOPO2v2c_i2_LSB.zip'

## ETOPO5
# see: https://www.ngdc.noaa.gov/mgg/global/etopo5.HTML
#
#ETOPO5 was generated from a digital data base of land and sea-floor elevations on a 5-minute latitude/longitude grid. The resolution of the gridded data varies from true 5-minute for the ocean floors, the USA., Europe, Japan,and Australia to 1 degree in data-deficient parts of Asia, South America, northern Canada, and Africa.
#
# note: not zip format, direct data file - case not handled yet, will used downsampling from etopo2
# (deprecated)
#filename_web5 = 'ETOPO5.DAT'
#url_etopo5 = 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO5/TOPO/ETOPO5/ETOPO5.DAT'

## MARS
# see: https://pds-geosciences.wustl.edu/missions/mgs/megdr.html
#
# MOLA topography data set
# 32-pixel per degree = 0.03125 degree resolution ~ 1.875 arc-minute
filename_webMOLA2lbl = 'megt90n000fb.lbl'
filename_webMOLA2 = 'megt90n000fb.img'
url_MOLA2 = 'https://pds-geosciences.wustl.edu/mgs/mgs-m-mola-5-megdr-l3-v1/mgsl_300x/meg032/megt90n000fb.img'


# plots PNM image
plot_image = True

####################################################################


def check_status(status):
    if status != 0:
        print("error: status returned ",status)
        sys.exit(status)
    return


def get_url(topo):
    """
    returns url for data download
    """
    global url_etopo1,url_etopo2 #url_etopo5
    # sets url
    url = ''
    if topo == 'etopo1':
        url = url_etopo1
    elif topo == 'etopo2':
        url = url_etopo2
    elif topo == 'etopo4':
        url = url_etopo2 # will get re-sampled
    elif topo == 'etopo5':
        url = url_etopo2 # will get re-sampled
    elif topo == 'etopo15':
        url = url_etopo2 # will get re-sampled
    elif topo == 'mars1.875':
        url = url_MOLA2
    elif topo == 'mars2':
        url = url_MOLA2  # will get re-sampled
    elif topo == 'mars4':
        url = url_MOLA2  # will get re-sampled
    else:
        print("Invalid topo argument %s, only recognizes etopo1,etopo2,etopo4,etopo5,etopo15,mars2,mars4" % topo)
        sys.exit(1)
    return url


def download_data_file(topo,filename):
    """
    downloads topo/bathy file from web
    """
    global filename_webMOLA2,filename_webMOLA2lbl

    print("download:")

    # checks if file already downloaded
    if os.path.isfile(filename):
        print("  file %s exists already, skipping download" % filename)
        return

    # download from web
    url = get_url(topo)
    print("url: %s" % url)

    # downloads large files as stream
    r = requests.get(url,stream=True)

    # file status
    if r.status_code == 200:
        status = 'success'
    elif r.status_code == 404:
        status = 'file not found'
    else:
        status = str(r.status_code)

    # file length
    length = r.headers.get('content-length')
    if length is None: length = 0

    # info
    print("url file info:")
    print("  status  : ",status)
    print("  header  : ",r.headers['content-type'])
    print("  encoding: ",r.encoding)
    print("  length  : {} {}".format(int(length) / 1024.0 / 1024.0,"(MB)"))

    # check
    if r.status_code == 404:
        print("Error: file not found",r.status_code)
        sys.exit(1)

    # downloads data file
    if 'mars' in topo:
        # data file in .img format
        filename = filename_webMOLA2
    else:
        # data file downloaded as .zip file
        filename = 'topo.zip'

    # saves data file as topo.zip
    with open(filename, 'wb') as f:
        if length == 0:
            # no content length header
            f.write(r.content)
        else:
            # status bar
            print("start downloading...")
            incr = 0
            t_start = timeit.default_timer()
            for data in r.iter_content(chunk_size=4096):
                incr += len(data)
                f.write(data)
                done = int(100.0 * incr / int(length))
                # download speed MB/s
                t_elapsed = timeit.default_timer() - t_start
                total_data = incr / 1024. / 1024.
                bandwidth = total_data / t_elapsed
                # status bar
                sys.stdout.write("\r[%s%s] %6.1f MB  %4.2f MB/s" % ('=' * done, ' ' * (100-done), total_data, bandwidth ) )
                sys.stdout.flush()
        print("")
        print("written: ",f.name)

    # extract data file
    if 'mars' in topo:
        # already data format
        # downloads also table info file (descriptor to actual data file, needed for later file conversions with GDAL/GMT)
        lbl_file = requests.get(filename_webMOLA2lbl)
        open(filename_webMOLA2,'w').write(lbl_file.content)
    else:
        # unzip
        with zipfile.ZipFile(filename, 'r') as zipObj:
            # list
            print("zip files: ",zipObj.namelist())
            # Extract all the contents of zip file in current directory
            print("zip extracting file...")
            zipObj.extractall()
    print("")

def resample_topo_file(topo,filename_web,filename):
    """
    re-samples etopo2 to a coarser grid format
    """
    # re-sampling
    print("resampling:")

    # sampling interval
    if topo == 'etopo1':
        Ival = 1   # to 1 minute resampling (resampling as .tif can lead to problems)
    elif topo == 'etopo2':
        Ival = 2   # to 2 minute resampling
    elif topo == 'etopo4':
        Ival = 4   # to 4 minute downsampling
    elif topo == 'etopo5':
        Ival = 5   # to 5 minute downsampling
    elif topo == 'etopo15':
        Ival = 15  # to 15 minute downsampling
    elif topo == 'mars1.875':
        Ival = 1.875
    elif topo == 'mars2':
        Ival = 2  # to 2 minute resampling
    elif topo == 'mars4':
        Ival = 4  # to 2 minute resampling
    else:
        print("invalid topo parameters for resampling, exiting...")
        sys.exit(1)

    print("  resampling to {}m".format(Ival))
    if 'mars' in topo:
        gridfile = "MARSTOPO{}.grd".format(Ival)
    else:
        gridfile = "ETOPO{}.grd".format(Ival)

    if gridfile != filename:
        print("inconsistent filenames: gridfile {} and filename {}".format(gridfile,filename))
        print("please check script, exiting...")
        sys.exit(1)

    # checks if file already done
    if os.path.isfile(filename):
        print("  file %s exists already, skipping resample" % filename)
    else:
        # sampling interval
        interval = "-I{}m".format(Ival)
        grid = "-G{}=ns".format(filename)

        if topo == 'etopo1':
            # converting from raw binary to grid file
            cmd = 'xyz2grd '   + filename_web + ' -V -Rd ' + grid + ' -r -ZTLh ' + interval

        elif topo == 'etopo2' or topo == 'etopo4' or topo == 'etopo5' or topo == 'etopo15':
            if filename_web == 'ETOPO2v2c_i2_LSB.bin':
                # converting from etopo2v2c raw binary to 2m grid file first
                print("  converting binary file to grid file: etopo2v2c.grd")
                cmd = 'xyz2grd ETOPO2v2c_i2_LSB.bin -Rd -I2m -Getopo2v2c.grd -r -ZTLh -V'
                print("> ",cmd)
                status = subprocess.call(cmd, shell=True)
                check_status(status)
                # resample
                cmd = 'grdsample etopo2v2c.grd -V -Rd ' + grid + ' ' + interval
            else:
                # direct resample
                cmd = 'grdsample ' + filename_web + ' -V -Rd ' + grid + ' ' + interval

        elif 'mars' in topo:
            # converting to raw binary grid file
            # converts to GeoTiff
            # > gdal_translate -of GTiff megt90n000fb.lbl topo.tif
            # converts to GeoTiff, method 3
            # see: https://pds-geosciences.wustl.edu/mgs/mgs-m-mola-5-megdr-l3-v1/mgsl_300x/extras/mola_envi.pdf
            # > gdal_translate -a_ullr upperLeftX UpperLeftY LowerRightX LowerLeftY megt90n000fb.lbl topo.tif
            # converts to ETOPO file
            # > xyz2grd megt90n000fb.img -Gtopo.grd=ns -R0/360/-90/90 -I1.875m -ZTLhw -fg -r -Vl
            # plot image
            # > grdimage topo.grd -JX6i+ -I+d -P -Vl > topo.ps
            print("  converting binary file to grid file: "+filename_web)
            cmd = 'xyz2grd '   + filename_web + ' -Gmarstopo1.875.grd -R0/360/-90/90 -I1.875m -r -ZTLhw -fic -fog -V'
            print("> ",cmd)
            status = subprocess.call(cmd, shell=True)
            check_status(status)
            # resample
            cmd = 'grdsample marstopo1.875.grd -V -Rd ' + grid + ' ' + interval
        else:
            print("invalid topo %s for resampling" % topo)
            sys.exit(1)

        print("  creating grid file: {}".format(gridfile))
        print("> ",cmd)
        status = subprocess.call(cmd, shell=True)
        check_status(status)

    print("")


def extract_to_ascii(topo,filename_in,filename_out):
    """
    takes extracted binary file and outputs topo values in ascii-format
    """
    print("extracting elevation:")

    # checks if file already downloaded
    if os.path.isfile(filename_out):
        print("  file %s exists already, skipping extraction" % filename_out)
    else:
        # file info
        cmd = 'grdinfo ' + filename_in
        print("> ",cmd)
        status = subprocess.call(cmd, shell=True)
        check_status(status)
        print("")

        # converts to ascii, extracts elevation values
        print("  extracting elevation values to ascii file...")
        cmd = 'grd2xyz -Rg ' + filename_in + ' | awk \'{ print $3 }\' > ' + filename_out
        print("> ",cmd)
        status = subprocess.call(cmd, shell=True)
        check_status(status)
        print("written to: %s" % filename_out)

    print("")
    print("  gmt info:")
    cmd = 'gmtinfo ' + filename_out
    print("> ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")


def create_padded_array(arr,n=3):
    """
    padding a 2D-array arr, uses periodic boundary along x-direction, edge boundary along y
    """
    # x-axis padding with wrap values (re-starting e/w)
    a = np.pad(arr,((n,n),(0,0)),mode='wrap')
    # y-axis padding with edge values (north / south poles)
    arr_padded = np.pad(a,((0,0),(n,n)),mode='edge')
    #arr = np.roll(np.roll(arr,shift=-x+1,axis=0),shift=-y+1,axis=1)
    return arr_padded


def sum_neighbors_with_padded_array(arr_padded,x,y,n=3):
    """
    returns the sum of the neighbors
    """
    # assumes arr is a padded array with padding in x/y direction of size n
    # shift index to have center at x,y
    ix = x + n
    iy = y + n
    #debug
    #print(b)
    #for ix_b in np.arange(ix-n,ix+n+1):
    #    for iy_b in np.arange(iy-n,iy+n+1):
    #        print("array index ",ix_b,iy_b," topo ",ix,iy,b[ix_b,iy_b])
    #print(b[ix-n:ix+n+1,iy-n:iy+n+1])
    # sum to average
    total_sum = arr_padded[ix-n:ix+n+1,iy-n:iy+n+1].sum()
    return total_sum

def smooth_topo_bathy_fortran_tool(topo,filename_in,filename_out,SIZE_FILTER_ONE_SIDE):
    """
    smooths elevation by fortran tool
    """
    # checks if file already downloaded
    if os.path.isfile(filename_out):
        print("  file %s exists already " % filename_out)
    else:
        print("  smoothing by fortran tool")
        cmd = 'gfortran smooth_topo_bathy_PPM_image.f90'
        print("> ",cmd)
        status = subprocess.call(cmd, shell=True)
        check_status(status)
        print("")

        cmd = './a.out {} {}'.format(filename_in,SIZE_FILTER_ONE_SIDE)
        print("> ",cmd)
        status = subprocess.call(cmd, shell=True)
        check_status(status)
        print("")

def smooth_topo_bathy(topo,filename_in,filename_out,SIZE_FILTER_ONE_SIDE=1,NX_BATHY=0,NY_BATHY=0):
    """
    smooths elevation by taking average value within neighboring data points
    """
    print("smoothing:")
    # pre-defined dimensions
    print("  topo/bathy dimensions:")
    print("  NX_BATHY = ",NX_BATHY)
    print("  NY_BATHY = ",NY_BATHY)
    print("")

    # checks if file already downloaded
    if os.path.isfile(filename_out):
        print("  file %s exists already, reading in..." % filename_out)
        print("  (this might take a while)...")
        ibathy_topo = np.loadtxt(filename_out,dtype='int')
        ibathy_topo = np.reshape(ibathy_topo, (NX_BATHY,NY_BATHY), order='F')
        return ibathy_topo

    # reads in data from ascii file
    print("  loading data from file: %s" % filename_in)
    print("  (this might take a while)...")
    data = np.loadtxt(filename_in)
    print("")

    # checks number of entries
    total_size = data.size
    print("  number of total data points = ",total_size)
    if total_size != NX_BATHY * NY_BATHY:
        print("Error dimension: size for %s should be %d" % (topo,NX_BATHY*NY_BATHY))
        sys.exit(1)

    # fortran shaped array as output format from gmt command grd2xyz
    topo_bathy = np.reshape(data, (NX_BATHY,NY_BATHY), order='F')
    # or explicit reading
    #topo_bathy = np.zeros((NX_BATHY,NY_BATHY))
    # fortran (column-major) order
    #for iy in np.arange(NY_BATHY):
    #    for ix in np.arange(NX_BATHY):
    #        # E/W-direction corresponds to x
    #        # N/S-direction             to y
    #        topo_bathy[ix,iy] = data[ix+iy*NX_BATHY]

    #debug
    #NX_BATHY = 4
    #NY_BATHY = 2
    #topo_bathy = np.arange(4 * 2).reshape((4, 2))
    #print(topo_bathy)

    print("")
    print("min and max of topography before smoothing = %f / %f" %(topo_bathy.min(),topo_bathy.max()))
    print("")

    # smoothed array
    ibathy_topo = np.zeros((NX_BATHY,NY_BATHY),dtype='int')

    # creates ibathy_topo array
    if SIZE_FILTER_ONE_SIDE == 0:
        print("  no smoothing...")
        # explicit looping
        for ix_current in np.arange(NX_BATHY):
            for iy_current in np.arange(NY_BATHY):
                ibathy_topo[ix_current,iy_current] = int(topo_bathy[ix_current,iy_current])


    else:
        # smoothing: average value in x/y neighborhood

        # see also: DATA/topo_bathy/smooth_topo_bathy_PPM_image.f90
        print("  using SIZE_FILTER_ONE_SIDE: {}".format(SIZE_FILTER_ONE_SIDE))
        # checks filter size
        #if SIZE_FILTER_ONE_SIDE < 1:
        #    print("SIZE_FILTER_ONE_SIDE must be greater than 1 for filter")
        #    sys.exit(1)
        print("  size of window filter is (%d x %d)" %(2*SIZE_FILTER_ONE_SIDE+1,2*SIZE_FILTER_ONE_SIDE+1))
        print("")
        area_window = np.float((2*SIZE_FILTER_ONE_SIDE+1)**2)


        count = 0
        incr = int(total_size / 100.0)

        print("  computing smoothed values...")
        t_start = timeit.default_timer()
        if 0 == 1:
            # explicit looping
            for ix_current in np.arange(NX_BATHY):
                for iy_current in np.arange(NY_BATHY):
                    # loop on points in window to compute sum
                    # compute min and max of window
                    ix_min = ix_current - SIZE_FILTER_ONE_SIDE
                    ix_max = ix_current + SIZE_FILTER_ONE_SIDE

                    iy_min = iy_current - SIZE_FILTER_ONE_SIDE
                    iy_max = iy_current + SIZE_FILTER_ONE_SIDE

                    value_sum = 0.0
                    #debug
                    #b = np.zeros((2*SIZE_FILTER_ONE_SIDE+1,2*SIZE_FILTER_ONE_SIDE+1))

                    for ix in np.arange(ix_min,ix_max+1):
                        for iy in np.arange(iy_min,iy_max+1):
                            # copy current value
                            ix_value = ix
                            iy_value = iy

                            # avoid edge effects, use periodic boundary in Xmin and Xmax
                            if ix_value < 0: ix_value = ix_value + NX_BATHY
                            if ix_value > NX_BATHY-1: ix_value = ix_value - NX_BATHY

                            # avoid edge effects, use rigid boundary in Ymin and Ymax
                            # *not* periodic, because South and North poles must not be merged
                            if iy_value < 0: iy_value = 0
                            if iy_value > NY_BATHY-1: iy_value = NY_BATHY-1

                            # compute sum
                            value_sum = value_sum + topo_bathy[ix_value,iy_value]

                            #debug
                            #ix_b = ix-ix_current+SIZE_FILTER_ONE_SIDE
                            #iy_b = iy-iy_current+SIZE_FILTER_ONE_SIDE
                            #print("array index ",ix_b,iy_b," topo ",ix,iy,topo_bathy[ix_value,iy_value])
                            #b[ix_b,iy_b] = topo_bathy[ix_value,iy_value]

                    #debug
                    #print(b)
                    #print("values:",ix_current,iy_current,value_sum)
                    #sys.exit()

                    # assign mean value to filtered point
                    ibathy_topo[ix_current,iy_current] = int(np.around(value_sum / area_window))

                    # progress info output
                    count += 1
                    if count % incr == 0:
                        done = int(100.0 * count/total_size)  # in percent
                        # timing
                        t_elapsed = timeit.default_timer() - t_start
                        # status bar
                        #sys.stdout.write("\r[%s%s]" % ('=' * done, ' ' * (100-done)) )
                        sys.stdout.write("\r{:3d} % done / elapsed time {:8.1f} s / estimated total time {:8.1f} s".format(done,t_elapsed,t_elapsed/done*100.0))
                        sys.stdout.flush()

        else:
            # using padded array
            # padding edges
            print("  padding array for corresponding boundary effects in E/W and N/S direction")
            topo_bathy = create_padded_array(topo_bathy,n=SIZE_FILTER_ONE_SIDE)

            for ix_current in np.arange(NX_BATHY):
                for iy_current in np.arange(NY_BATHY):
                    # gets sum of temporary neighbors array around x,y position
                    value_sum = sum_neighbors_with_padded_array(topo_bathy,ix_current,iy_current,n=SIZE_FILTER_ONE_SIDE)

                    #debug
                    #print("values:",ix_current,iy_current,value_sum)
                    #sys.exit()

                    # assign mean value to filtered point
                    ibathy_topo[ix_current,iy_current] = int(np.around(value_sum / area_window))

                    # progress info output
                    count += 1
                    if count % incr == 0:
                        done = int(100.0 * count/total_size)  # in percent
                        # timing
                        t_elapsed = timeit.default_timer() - t_start
                        # status bar
                        #sys.stdout.write("\r[%s%s]" % ('=' * done, ' ' * (100-done)) )
                        sys.stdout.write("\r{:3d} % done / elapsed time {:8.1f} s / estimated total time {:8.1f} s".format(done,t_elapsed,t_elapsed/done*100.0))
                        sys.stdout.flush()

    print("")

    min_val = ibathy_topo.min()
    max_val = ibathy_topo.max()
    print("")
    print("min and max of topography after smoothing = %f / %f" %(min_val,max_val))
    print("")

    print("  saving the smoothed model: %s" % filename_out)
    with open(filename_out,'w') as f:
        # saves in fortran (column-major) order
        for iy in np.arange(NY_BATHY):
            for ix in np.arange(NX_BATHY):
                f.write("{}\n".format(ibathy_topo[ix,iy]))

    return ibathy_topo

def plot_PPM_image(ibathy_topo,filename='image_topo_bathy.pnm'):
    """
    plots a (gray) image of the given topo/bathy array
    """
    # PPM image: create image with grey levels
    print("creating PNM image:")

    image_data = np.copy(ibathy_topo)

    (xdim,ydim) = np.shape(image_data)
    min_val = image_data.min()
    max_val = image_data.max()

    print("  image dimensions: {} x {}".format(xdim,ydim))
    print("  topography min/max values = %f / %f" %(min_val,max_val))

    # scales image values between [0,255]
    for ix in np.arange(xdim):
        for iy in np.arange(ydim):
            image_data[ix,iy] = 255 * (image_data[ix,iy] - min_val) / (max_val - min_val)
            if image_data[ix,iy] < 1: image_data[ix,iy] = 1
            if image_data[ix,iy] > 255: image_data[ix,iy] = 255

    # store image in PNM format with grey levels
    # create the PNM image
    with open(filename,'w') as f:
        # creating the header
        f.write('P3\n')
        f.write('{:6d} {:6d}\n'.format(xdim,ydim))
        f.write('255\n')
        # image values
        if 0 == 1:
            # writes values (red == green == blue) to produce grey levels
            # saves in fortran (column-major) order
            for iy in np.arange(ydim):
                for ix in np.arange(xdim):
                    val = image_data[ix,iy]
                    f.write('{:3d} {:3d} {:3d}\n'.format(val,val,val))
        else:
            # fun coloring
            # saves in fortran (column-major) order
            for iy in np.arange(ydim):
                for ix in np.arange(xdim):
                    val = image_data[ix,iy]
                    # actual elevation
                    elevation = val / 255.0 * (max_val - min_val) + min_val
                    # in [0,1]
                    c = val / 255.0
                    # power-scaling for nicer effect
                    c = c*c
                    color = int(c * 255.0)
                    # coloring
                    if elevation < 0.0:
                        # blueish - below sea-level
                        f.write('{:3d} {:3d} {:3d}\n'.format(int(val*0.5),int(val*0.5),val))
                    elif elevation <= 1500.0:
                        # greenish - vegetation
                        f.write('{:3d} {:3d} {:3d}\n'.format(int(color*0.5),color,int(color*0.2)))
                    else:
                        # gray-white - snow in high altitudes
                        f.write('{:3d} {:3d} {:3d}\n'.format(color,color,color))

        print("  image written to: image_topo_bathy.pnm")
    print("")


def convert_data_to_specfem_format(ibathy_topo,filename_out):
    """
    creates binary file for SPECFEM mesher
    """
    print("converting to specfem format:")

    # re-shape 2d-array ibathy_topo to a 1-d array
    # (fortran shaped array as output format)
    (xdim,ydim) = np.shape(ibathy_topo)
    data = np.reshape(ibathy_topo, (xdim*ydim), order='F')
    # or explicit
    #data = np.zeros(xdim*ydim)
    # fortran (column-major) order
    #for iy in np.arange(ydim):
    #    for ix in np.arange(xdim):
    #        # E/W-direction corresponds to x
    #        # N/S-direction             to y
    #        data[ix+iy*xdim] = ibathy_topo[ix,iy]

    # Convert to 16-bit integers
    data2 = data.astype(np.int16)
    if any(data != data2):
        print('Warning: Data set does not fit in signed 16-bit integers!')

    # Output file
    with open(filename_out, 'wb') as outf:
        # Add a byte-order mark
        byteorder = np.array([0x1234], dtype=np.int16)
        byteorder.tofile(outf)

        # Save output file
        data2.tofile(outf)

        print("  written to: %s" % filename_out)
        print("")

    # in case for transfering to other systems, one could further compress/decompress the binary file:
    # > bzip2 <topo-file>.bin
    # will create a file like <topo-file.bin.bz2 which can be decompressed with
    # > bunzip2 <topo-file.bin.bz2


def create_topo_bathy(topo):
    """
    downloads and creates SPECFEM3D_GLOBE readable topo file
    """
    global filename_web1,filename_web2,filename_webMOLA2,plot_image

    print("topo/bathy:")
    print("  etopo: %s" % topo)

    # sets filenames
    if topo == 'etopo1':
        INTERVAL = 1
        SIZE_FILTER_ONE_SIDE = 3
        NX_BATHY = 21600
        NY_BATHY = 10800
        filename_web = filename_web1

    elif topo == 'etopo2':
        INTERVAL = 2
        SIZE_FILTER_ONE_SIDE = 3
        NX_BATHY = 10800
        NY_BATHY = 5400
        filename_web = filename_web2

    elif topo == 'etopo4':
        INTERVAL = 4
        SIZE_FILTER_ONE_SIDE = 7
        NX_BATHY = 5400
        NY_BATHY = 2700
        filename_web = filename_web2

    elif topo == 'etopo5':
        INTERVAL = 5
        SIZE_FILTER_ONE_SIDE = 7
        NX_BATHY = 4320
        NY_BATHY = 2160
        filename_web = filename_web2

    elif topo == 'etopo15':
        INTERVAL = 15
        SIZE_FILTER_ONE_SIDE = 3
        NX_BATHY = 1440
        NY_BATHY = 720
        filename_web = filename_web2

    elif topo == 'mars1.875':
        INTERVAL = 1.875          # original 32-pixel per degree = 0.03125 degree -> 1.875 arc-minutes
        SIZE_FILTER_ONE_SIDE = 0  # no filtering
        NX_BATHY = 11520
        NY_BATHY = 5760
        filename_web = filename_webMOLA2

    elif topo == 'mars2':
        INTERVAL = 2              # original 32-pixel per degree = 0.03125 degree -> 1.875 arc-minutes, resample to 2m
        SIZE_FILTER_ONE_SIDE = 0  # no filtering
        NX_BATHY = 10800
        NY_BATHY = 5400
        filename_web = filename_webMOLA2

    elif topo == 'mars4':
        INTERVAL = 4              # resample to 2m
        SIZE_FILTER_ONE_SIDE = 3
        NX_BATHY = 5400
        NY_BATHY = 2700
        filename_web = filename_webMOLA2

    else:
        print("Unrecognized option %s" % topo)
        sys.exit(1)

    # 1-degree in km
    if 'mars' in topo:
        degree_distance = 2.0 * 3390.0 * 3.1415926 / 360.0
    else:
        degree_distance = 2.0 * 6371.0 * 3.1415926 / 360.0
    # distance for minutes interval
    grid_resolution_minutes = INTERVAL * degree_distance / 60.0

    print("  interval: {} minutes".format(INTERVAL))
    print("  interval grid point spacing: {:6.2f} km".format(grid_resolution_minutes))
    print("  filter size one side: {}".format(SIZE_FILTER_ONE_SIDE))
    print("  filter size one side distance: {:6.2f} km".format(SIZE_FILTER_ONE_SIDE*grid_resolution_minutes))
    print("  topo_bathy size: ({} x {})".format(NX_BATHY,NY_BATHY))
    print("")

    # GMT grid filename
    if 'mars' in topo:
        filename_grid = "MARSTOPO{}.grd".format(INTERVAL)
        # ASCII filenames
        filename_out = "topo_bathy_marstopo{}_unmodified_unsmoothed.dat".format(INTERVAL)
        filename_out_smoothed = "topo_bathy_marstopo{}_smoothed_window_{}.dat".format(INTERVAL,SIZE_FILTER_ONE_SIDE)
    else:
        filename_grid = "ETOPO{}.grd".format(INTERVAL)
        # ASCII filenames
        filename_out = "topo_bathy_etopo{}_unmodified_unsmoothed.dat".format(INTERVAL)
        filename_out_smoothed = "topo_bathy_etopo{}_smoothed_window_{}.dat".format(INTERVAL,SIZE_FILTER_ONE_SIDE)
        #filename_grid = 'etopo1_ice_c_resampled_at_2minutes.grd'
        #filename_out = 'topo_bathy_etopo1_ice_c_resampled_at_2minutes_original_unmodified_unsmoothed.dat'

    # data file download
    download_data_file(topo,filename_web)

    # re-samples topo
    resample_topo_file(topo,filename_web,filename_grid)

    # extracts elevation
    extract_to_ascii(topo,filename_grid,filename_out)

    # smoothing
    if topo == 'etopo1' or topo == 'etopo2':
        # etopo1/etopo2 smoothing
        # takes too long in python, using provided fortran tool
        smooth_topo_bathy_fortran_tool(topo,filename_out,filename_out_smoothed,SIZE_FILTER_ONE_SIDE)

        # conversion to SPECFEM binary
        # taken from file: convert_etopo_files_from_specfem_ASCII_to_binary_topo_format_ascii2bin.py
        # uses conversion by chunks, better suited for large files
        from convert_etopo_files_from_specfem_ASCII_to_binary_topo_format_ascii2bin import convert_etopo_ascii2bin
        convert_etopo_ascii2bin(filename_out_smoothed,filename_out_smoothed + '.bin')

    else:
        # uses python routine for smoothing
        ibathy_topo = smooth_topo_bathy(topo,filename_out,filename_out_smoothed,SIZE_FILTER_ONE_SIDE,NX_BATHY,NY_BATHY)

        # PPM image
        if plot_image: plot_PPM_image(ibathy_topo)

        # conversion to SPECFEM binary
        convert_data_to_specfem_format(ibathy_topo,filename_out_smoothed + '.bin')

    print("all done")


if __name__ == '__main__':
    # gets arguments
    if '--help' in sys.argv or '-h' in sys.argv or len(sys.argv) != 2:
        print("usage: ./create_topo_bathy_file.py topo [== etopo1,etopo2,etopo4,etopo5,etopo15,mars2,mars4]")
        sys.exit(1)
    else:
        topo = sys.argv[1]

    create_topo_bathy(topo)

