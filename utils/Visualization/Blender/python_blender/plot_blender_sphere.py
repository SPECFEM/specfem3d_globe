#!/usr/bin/env blender --background --python-use-system-env --python plot_blender_sphere.py --
#
# or for example on Mac:
#  /usr/bin/env /Applications/Blender.app/Contents/MacOS/Blender --background --python-use-system-env --python plot_blender_sphere.py --
#
#
# run with: > ./plot_blender_sphere.py --help
#
###############################################################################################

import sys
import os
import glob
import time


# import user scripts
from moviedata import moviedata
from blender_sphere import blender_sphere

###############################################################################################
## USER parameters

## renderer
# render image size
blender_img_resolution_X = 1600
blender_img_resolution_Y = 1600

# blender engine
blender_engine = 'BLENDER_EEVEE'  # 'CYCLES', 'BLENDER_EEVEE'
blender_device = 'CPU'            # for cycles engine: 'CPU', 'GPU', ..

# suppressing renderer output
suppress_renderer_output = False

## textures
# surface image
earth_image = "resources/world.jpg"     # "resources/earth.jpg"
# at night
earth_night = "resources/night.jpg"
# clouds
earth_clouds = "resources/clouds.jpg"
# topography
earth_topo  = "resources/topo_8192.png"

## moviedata
# moviedata texture grid resolution
moviedata_grid_resolution = 2 * max(blender_img_resolution_X,blender_img_resolution_Y)

# color power scaling (0.0==turn off)
moviedata_color_power_scaling = 0.3

# fixes colormap maximum value
moviedata_colormap_max = None

# smoothing
moviedata_smooth = False

## animation
# camera position field of view (fov)
field_of_view = 40.0

# animation rotation rate (in degrees)
animation_rotation_rate = - 0.2

###############################################################################################


# main routine
def plot_blender_sphere(fov=50.0,animation=False,datafile=""):
    """
    renders image for (earth) sphere with textures
    """
    # set current directory, in case we need it to load files
    dir = os.getcwd()
    print("current directory: ",dir)
    print("")

    # texture files
    filename_image = ""
    filename_topo = ""
    filename_night = ""
    filename_clouds = ""
    if earth_image : filename_image  = dir + "/" + earth_image
    if earth_topo  : filename_topo   = dir + "/" + earth_topo
    if earth_night : filename_night  = dir + "/" + earth_night
    if earth_clouds: filename_clouds = dir + "/" + earth_clouds

    # gets number of input moviedata files (wildcard possible: moviedata00*)
    datafiles = []
    num_datafiles = 0
    if datafile:
        print("moviedata file: ",datafile)
        # check if multiple files specified with wildcard (e.g., moviedata00*)
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

        # number of moviedata files
        num_datafiles = len(datafiles)
        print("  number of moviedata files = ",num_datafiles)

        # checks if file(s) exists
        for filename in datafiles:
            # checks if file exists
            if not os.path.isfile(filename):
                print("Please check if file exists: ",filename)
                sys.exit(1)
        print("  files ok")
        print("")

    # setup blender renderer scene
    print("creating blender scene... ")
    renderer = blender_sphere(verbose=True,
                              texture_globe=filename_image,
                              texture_clouds=filename_clouds,
                              texture_night=filename_night,
                              texture_topo=filename_topo,
                              animation=animation,
                              animation_rotation_degree=animation_rotation_rate)

    # sets rendering defaults
    renderer.blender_engine = blender_engine
    renderer.blender_device = blender_device
    renderer.blender_img_resolution_X = blender_img_resolution_X
    renderer.blender_img_resolution_Y = blender_img_resolution_Y

    # adds scene effects
    renderer.add_scene_effects()

    # adds lights
    renderer.add_sun()

    # adds moon lights
    renderer.add_moon()

    # adds camera
    renderer.add_camera()

    # adds globe sphere
    renderer.add_sphere()

    # adds clouds layer
    renderer.add_clouds()

    # read/import moviedata and create texture shader
    if num_datafiles > 0:
        # moviedata helper
        movdata = moviedata(verbose=True)

        # output data norm-component (1==Z,2==N,3==E,4==norm)
        movdata.use_component = 4

        # smoothing
        if moviedata_smooth:
            movdata.use_smoothing = moviedata_smooth
            movdata.use_smoothing_kernel_size = 3
            if moviedata_grid_resolution > 1000:
                movdata.use_smoothing_kernel_size = 10

        # reads in moviedata and creates texture images
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

            # timing
            tic = time.perf_counter()

            # for example: moviedata000200
            lat,lon,data = movdata.read_moviedata(filename,verbose=True)

            # matplotlib/scipy linear interpolation
            xi,yi,data_gridded = movdata.interpolate_data(lat,lon,data,
                                                          grid_resolution=moviedata_grid_resolution)

            # creates texture image from moviedata
            renderer.create_moviedata_textureImage(data_gridded,movdata,i,
                                                   power_scaling=moviedata_color_power_scaling,
                                                   colormap_max=moviedata_colormap_max)

            # timing
            toc = time.perf_counter()
            print("elapsed time for creating moviedata texture is {:0.4f} seconds\n".format(toc - tic))

        # adds moviedata image as texture node
        renderer.add_moviedata(movdata)

    # rendering
    print("***************************************************")
    print("")
    print("rendering")
    print("")
    print("***************************************************")
    print("")

    # rendering
    if animation:
        # setup animation keyframes
        renderer.setup_animation_frames(num_datafiles)

        # renders animation
        renderer.output_animation(dir,fov,suppress=suppress_renderer_output)

    else:
        # still image rendering
        if datafile:
            # loop over moviedata files
            i = 0
            for filename in datafiles:
                i += 1
                # sets moviedata frames
                renderer.setup_image_frame(i,num_datafiles)

                # name appendix
                name_appendix = os.path.basename(filename)

                # renders image
                renderer.output_image(dir,fov,appendix=name_appendix,suppress=suppress_renderer_output)
        else:
            # renders image
            renderer.output_image(dir,fov,suppress=suppress_renderer_output)

    # all done
    print("rendering done")
    print("")


def usage():
    print("usage: ./plot_blender_sphere.py [--anim] [--anim-rotation-rate=val] [--color-max=val]")
    print("                                [--engine='type'] [--device='type'] [--fov=val]")
    print("                                [--moviedata=file] [--moviedata-color-power-scaling=val] [--moviedata-smooth]")
    print("                                [--no-clouds] [--suppress] ")
    print("  with")
    print("     --anim                  - renders animation, rotates sphere")
    print("     --anim-rotation-rate    - sets animation rotation rate in degrees (default=-0.2)")
    print("     --color-max             - fixes maximum value of colormap for moviedata to val, e.g., 1.e-7)")
    print("     --device                - select render device 'CPU' or 'GPU' (default 'CPU')")
    print("     --engine                - select render engine 'CYCLES' or 'BLENDER_EEVEE' (default 'BLENDER_EEVEE')")
    print("     --fov                   - value for field of view (default=50)")
    print("     --moviedata             - file for moviedata, e.g., --moviedata=OUTPUT_FILES/moviedata003600")
    print("     --moviedata-color-power-scaling  - color power scaling factor for moviedata")
    print("     --moviedata-smooth      - turns on smoothing for moviedata")
    print("     --no-clouds             - turns off clouds (default turned on)")
    print("     --suppress              - turns off renderer stdout output (default turned on)")
    sys.exit(1)


if __name__ == '__main__':
    # init
    fov = field_of_view
    animation = False
    moviedatafile = ""

    # reads arguments
    #print("\nnumber of arguments: " + str(len(sys.argv)))
    i = 0
    for arg in sys.argv:
        i += 1
        #print("argument "+str(i)+": " + arg)
        # get arguments
        if "--help" in arg:
            usage()
        elif "--anim" in arg:
            animation = True
        elif "--anim-rotation-rate" in arg:
            animation = True
            animation_rotation_rate = float(arg.split('=')[1])
        elif "--color-max" in arg:
            moviedata_colormap_max = float(arg.split('=')[1])
        elif "--device" in arg:
            blender_device = arg.split('=')[1]
        elif "--engine" in arg:
            blender_engine = arg.split('=')[1]
        elif "--fov" in arg:
            fov = float(arg.split('=')[1])
        elif "--moviedata=" in arg:
            moviedatafile = arg.split('=')[1]
        elif "--moviedata-color-power-scaling" in arg:
            moviedata_color_power_scaling = float(arg.split('=')[1])
        elif "--moviedata-smooth" in arg:
            moviedata_smooth = True
        elif "--no-clouds" in arg:
            earth_clouds = ""
        elif "--suppress" in arg:
            suppress_renderer_output = True
        elif "--small" in arg:
            blender_img_resolution_X = 200
            blender_img_resolution_Y = 200
        elif i >= 8:
            print("argument not recognized: ",arg)

    # main routine
    plot_blender_sphere(fov,animation,moviedatafile)
