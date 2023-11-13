# Blender scripting example

Blender for 3D graphics creation<br>
https://www.blender.org/


## Installation

The python script `blender_sphere.py` uses Blender's python module. To use the script, we also need some interpolation routines from  matplotlib/scipy which are not provided in the default python version that blender internally uses. One possibility to use the systems python frameworks is to set an environment variable `BLENDER_SYSTEM_PATH`. For example, on Mac having python installed through [MacPorts](https://www.macports.org), one can set
```
export BLENDER_SYSTEM_PYTHON='/opt/local/Library/Frameworks/Python.framework/Versions/3.10/'
```
For this to work, the python version must match the internal python version from Blender. In this example, Blender version 3.6 uses a python version 3.10.

Another option is to install matplotlib into the provided Blender python version. For example, on Mac this can be used:
```
/Applications/Blender.app/Contents/Resources/3.6/python/bin/python3.10 -m pip install matplotlib
```


## Simulation setup

First, run a (global) simulation example, e.g., in `EXAMPLES/global_small/`, turn on the surface movie flag in the `DATA/Par_file` and set the output type to velocities:
   ```
   MOVIE_SURFACE                   = .true.
   ..
   MOVIE_VOLUME_TYPE               = 6
   ```

   This will create a series of `moviedata***` files in the example `OUTPUT_FILES/` folder.
   Create a (symbolic) link to `OUTPUT_FILES` into this rendering working directory, e.g.,
   ```
   ln -s ../../../../EXAMPLES/global_small/OUTPUT_FILES/
   ```


## Visualization

To plot an image with python's matplotlib module of such a `moviedata***` file, you can type:
```
./plot_moviedata.py --moviedata=OUTPUT_FILES/moviedata002800 --show
```
This script just uses the default python framework with matplotlib/scipy to interpolate and plot an image, without any need of blender's functionality. For more options on this script type `./plot_moviedata.py --help` . It is meant to check and quickly visualize the moviedata files.


### Blender textures

For the example here to work, the blender script `plot_blender_sphere.py` will need a few texture images for the rendering of the sphere (earth image, topography, clouds, night). These can be created and put into the `resources/` folder. These textures are not provided here, as the file sizes can become quite large for high-res images and would bloat up the github repository.

The `shakemovie` repository uses textures as well, you can use the script `convert_tga_to_jpg.sh` in the `resources/` folder to convert those `*.ppm` and `*.tga` files to high-res JPEG and PNG files, which are formats readable by blender.

You could run the script also with no textures by setting the texture paths to empty strings in the user parameter section in file `plot_blender_sphere.py`.


### Blender rendering

Once the texture files are created, you can set the corresponding file names as user parameters in the python script `plot_blender_sphere.py`, then type:
 ```
 ./plot_blender_sphere.py --moviedata=OUTPUT_FILES/moviedata002800
 ```

 This will create an image `out.moviedata002800.jpg`. For more script options, you can type `./plot_blender_sphere.py --help` .


Finally, to create an animation, you would call the script with multiple moviedata files as input. It will render all data and rotate the sphere to create a movie:
 ```
 ./plot_blender_sphere.py --moviedata=OUTPUT_FILES/moviedata00**
 ```

Please be patient, the rendering for many moviedata files can take a while. The movie file is stored as `out.anim.mp4`.


Feel free to modify and contribute any improvements to these python scripts - and have fun with blender for scientific visualization!
