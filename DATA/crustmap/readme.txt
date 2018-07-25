CRUSTMAPS Version 0.0.0.1

---------------------------------------------------------------------------------------------------
HOW TO USE THE INCLUDED CRUSTMAPS:

to use CRUSTMAPS append "_crustmaps" to your model name. for example:
MODEL    = s20rts_crustmaps



SPECFEM3D_GLOBE is set to use the crustmaps provided in the SPECFEM3D_GLOBE/DATA/crustmap/ folder.

-------------------------------------------------------------------------------------------------
WHAT ARE CRUSTMAPS?

the newly provided crustmaps in SPECFEM3D_GLOBE are a set of maps to provided a high resolution
crust model in SPECFEM3D_GLOBE on specified gridpoints in the crust.
Crustmaps can be used to integrate different high resolution crust models in lower resolution
models like crust2.0.

There are 4 types of crustmaps:

topography crustmaps (e.g. eucrustt3.cmap - eucrustt7.cmap)
s-velocity crustmaps (e.g. eucrusts3.cmap - eucrusts7.cmap)
p-velocity crustmaps (e.g. eucrustp3.cmap - eucrustp3.cmap)
density    crustmaps (e.g. eucrustr3.cmap - eucrustr7.cmap)

these are simple ascii data files and you can change them easily. The current crustmap resolution
is set to 1/4 degree.
The Values are ordered on a grid that runs from 90 degrees lat (0 degrees lng to 180 degrees lng ) to -90 degrees lat.
Every line in a crustmap files sets all values for one specific latitude:

(length of one line = 360 * 4 entries)
(total number of lines = 180 * 4)

The value for (lat = 89.875, long = 0.125) is stored in entry (line 1,column 1).
you get the value that is north west of your (lat, long) point by using:

line = int(0.5+((90.0-lat)*CRUSTMAP_RESOLUTION))
column = int0.5+(xlng*CRUSTMAP_RESOLUTION)

Values in between points are interpolated.
the topography file sets the thickness of every layer (positive number, depth under next layer)
the s velocity, p velocity and density maps set the respective values for every (lat,lng) point in every layer
The current Layer model follows crust2.0 (7 Layers). The ocean and ice layer of crust2.0 are not used,
since they can't be used by the code. The crystalline crust is (rather arbitrarily) seperated
into 3 layers, upper, middle and lower crust.

The provided maps combine The moho topography of crust2.0 and Eucrust with the velocity values from Crust2.0
(Layers not provided by Eucrust are taken from crust2.0 and stretched or squeezed according to the new
topography)

Specfem will stretch the Moho as regarding to the lowest layer.
The first two layers are regarded sediment layers (see option INCLUDE SEDIMENTS in constants.h)
but except from this option they are not treated differently.

----------------------------------------------------------------------------------------------------

HOW TO CREATE YOUR OWN CRUSTMODEL:

To create your own crustmodel you have to provide your own set of crustmaps. (you can either change the
filenames in model_crustmaps.f90 --> subroutine read_general_crustmap() or exchange the given files)

Make sure you provide the same Resolution as the default Crustmaps (1/4 degree) or change the Variable
CRUSTMAPS_RESOLUTION in constants.h.in. You can also change the Variable NLAYERS_CRUSTMAPS and try to use
a different number of crustmap layers (not tested).

If you have a set of new local crust layers, and want to integrate it in e.g. crust 2.0. It is easiest to
use a matlab/octave function like griddata to interpolate between all your given data points.

If your crust model is not complete (e.g. you have only 2 layers and not 3 or no velocity or density data)
you should try to combine your data with crust2.0 to provide a complete set of values.
You can also just change the different topography maps and the old velocity etc.. values are assigned to your
new points.

IMPORTANT JACOBIAN PROBLEMS:

if you experience Jacobian problems while running the mesher, you have to

I. smooth your model:
go to line 479 in the model_crustmaps.f90 and try to understand and change the following.

II.
if that doesn't work, you have to change your mesh resolution (NEX,...)
provide a less dramatic crustmap or turn off moho stretching.



