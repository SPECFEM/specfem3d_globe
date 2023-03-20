#!/usr/bin/env python
#
# plots a model defined as GLL or binary file
# e.g. OUTPUT_FILES/moviedata000200
#
# requires xyz file OUTPUT_FILES/moviedata_xyz.bin
#
from __future__ import print_function

import sys
import os
import glob

import numpy as np

# import user script
from moviedata import moviedata

###############################################################################################
# USER PARAMETERS

# interpolation: linear (by matplotlib/scipy) or nearest neighbor (numpy-only version)
use_linear_interpolation = True

# smoothing
use_smoothing = False

# grid image resolution
grid_resolution = 400 # 1200

# colormaps: "RdBu_r", "twilight", "seismic", "rainbow", "jet"
colormap = "jet"

# fixes colormap maximum value
colormap_max = None

# figure size
figure_size = (16,8)

###############################################################################################


def plot_model_image(xi,yi,data_gridded,name,verbose=False,show=False):
    """
    plots model or kernel image
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap

    # 2d array
    extent = [xi.min(),xi.max(),yi.min(),yi.max()]

    # limit size
    total_max = abs(data_gridded).max()

    if verbose:
        print("plot image: name = ",name)
        print("plot image: extent xmin/xmax/ymin/ymax = ",extent)
        print("plot image: data shape = ",np.shape(data_gridded))
        print("plot image: data max = ",total_max)

    # plotting
    fig = plt.figure(figsize=figure_size)
    plt.clf()
    plt.suptitle(name)

    if True:
        # cylindrical spherical projection
        m = Basemap(projection='cyl',resolution='c',
                    llcrnrlat=-90,urcrnrlat=90,
                    llcrnrlon=-180,urcrnrlon=180)
        #m.bluemarble()
        m.drawcoastlines()
        #m.shadedrelief(scale=0.5)
        origin = 'lower'
    else:
        # default for extent would plot on range [0,xmax,0,ymax]
        extent = None
        origin = 'lower'

    # 2d grid plot
    #Xi, Yi = np.meshgrid(xi, yi)
    #m.pcolormesh(Xi, Yi, data_gridded, latlon=True, cmap=colormap)

    # image interpolation
    # https://matplotlib.org/stable/gallery/images_contours_and_fields/interpolation_methods.html
    interp = 'none'  # 'none', 'nearest', 'bilinear', 'bicubic', 'spline16', ..

    # data image plot with explicit colormap maximum
    use_total_max_plot = True

    if use_total_max_plot:
        # determines a useful maximum value
        if total_max < 1.e-3:
            if total_max != 0.0:
                #total_max = 1.0 * 10**(int(np.log10(total_max))-1)  # example: 2.73e-11 limits to 1.e-11
                total_max = 1.0 * 10**(int(np.log10(total_max))-2)  # example: 2.73e-11 limits to 1.e-12
            else:
                total_max = 0.0

        # checks if fixing maximum value
        if colormap_max:
            total_max = colormap_max

        if verbose:
            if colormap_max:
                print("plot image: color scale max = ",total_max," (fixed)")
            else:
                print("plot image: color scale max = ",total_max)

        # sets min/max
        plt.imshow(data_gridded, extent=extent, origin=origin, interpolation=interp, cmap=colormap,
                   vmax=total_max, vmin=-total_max)

    else:
        # w/out value limits
        plt.imshow(data_gridded, extent=extent, origin=origin, interpolation=interp, cmap=colormap)

    plt.colorbar(location='bottom',extend='both',fraction=0.04,pad=0.04)  # draw colorbar

    # save as JPEG file
    filename = name + ".png"
    plt.savefig(filename)
    print("plotted as: ",filename)

    # show plot
    if show: plt.show()


##########################


def plot_moviedata_file(datafile,show=False):
    """
    reads binary file together with positions x/y/z and plots model as image
    """
    # user output
    print("input file: ",datafile)
    print("")

    # check if multiple files specified with wildcard (e.g., moviedata00*)
    datafiles = []
    if '*' in datafile:
        print("  files expanded:")
        files = glob.glob(datafile)
        files.sort()
        for filename in files:
            if "." in filename:
                # skip file name (e.g. moviedata002200.png)
                continue
            else:
                # keep file
                print("  ",filename)
                datafiles.append(filename)
        print("")
    else:
        datafiles.append(datafile)
    print("  number of moviedata files: ",len(datafiles))
    print("")
    # checks if file(s) exists
    for filename in datafiles:
        # checks if file exists
        if not os.path.isfile(filename):
            print("Please check if file exists: ",filename)
            sys.exit(1)

    # be verbose
    verbose = True

    # moviedata helper
    movdata = moviedata(verbose=verbose)

    # output data norm-component (1==Z,2==N,3==E,4==norm)
    movdata.use_component = 4

    # smoothing
    if use_smoothing:
        movdata.use_smoothing = use_smoothing
        # set kernel size
        movdata.use_smoothing_kernel_size = 4

    # reads in moviedata
    i = 0
    for filename in datafiles:
        # reads in model file
        i += 1
        print("")
        print("***************************************************")
        print("file {} out of {}".format(i,len(datafiles)))
        print("")
        print("moviedata: ",filename)
        print("***************************************************")
        print("")

        # for example: moviedata000200
        lat,lon,data = movdata.read_moviedata(filename,verbose)

        # gridding to regular 2D grid
        if use_linear_interpolation:
            # matplotlib/scipy linear interpolation
            xi,yi,data_gridded = movdata.interpolate_data(lat,lon,data,grid_resolution=grid_resolution)
        else:
            # manual nearest neighbor
            # (numpy-only version, no dependency on matplotlib or scipy)
            xi,yi,data_gridded = movdata.interpolate_data_NNglobal(lat,lon,data,grid_resolution=grid_resolution)

        # plot image
        print("plotting image...")
        plot_model_image(xi,yi,data_gridded,filename,show=show,verbose=verbose)

        print("")

    print("")
    print("done")


##########################


def usage():
    print("Usage: ./plot_moviedata.py --moviedata=filename [--color-max=val] [--colormap='jet'] [--resolution=val]")
    print("                           [--show] [--smooth]")
    print("  with")
    print("    --moviedata=filename - file for moviedata, e.g. OUTPUT_FILES/moviedata000200")
    print("                           (requires to have OUTPUT_FILES/moviedata_xyz.bin available as well)")
    print("    --color-max=val      - (optional) fixes maximum value of colormap to val (e.g., 1.e-7)")
    print("    --colormap='jet'     - (optional) sets colormap, e.g., 'jet', 'RdBu_r','seismic' (default is 'jet')")
    print("    --resolution=val     - (optional) sets grid resolution value")
    print("    --show               - (optional) show matplot image, otherwise will just be safed as .png file")
    print("    --smooth             - (optional) turns on smoothing for moviedata")
    sys.exit(1)


##########################

if __name__ == '__main__':
    # initializes
    show = False
    moviedatafile = ""

    # gets arguments
    if len(sys.argv) < 2:
        usage()

    # reads arguments
    i = 0
    for arg in sys.argv:
        i += 1
        #print("arg: ",arg)
        # get arguments
        if "--moviedata=" in arg:
            moviedatafile = arg.split('=')[1]
        elif "--color-max" in arg:
            colormap_max = float(arg.split('=')[1])
        elif "--colormap" in arg:
            colormap = arg.split('=')[1]
        elif "--show" in arg:
            show = True
        elif "--smooth" in arg:
            use_smoothing = True
        elif "--resolution" in arg:
            grid_resolution = int(arg.split('=')[1])
        elif i >= 1:
            print("argument not recognized: ",arg)

    # checks data file name
    if len(moviedatafile) < 1:
        usage()

    plot_moviedata_file(moviedatafile,show=show)
