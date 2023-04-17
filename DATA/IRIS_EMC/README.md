# IRIS Earth Model Collaboration (EMC)

Earth models provided by IRIS
http://ds.iris.edu/ds/products/emc-earthmodels/


## Model Download

To download the netCDF files, check out the following webpages:
* Alaska model:<br>
  https://ds.iris.edu/ds/products/emc-alaskajointinversion_rfvphhv-1berg2020/

* Southern California model:<br>
  http://ds.iris.edu/ds/products/emc-cvm_h_v15_1/


## Code compilation

The code must be configured with the `--with-netcdf` or `--with-emc` option, like:
```
./configure --with-netcdf NETCDF_FCFLAGS="<additional_fortran_flags> NETCDF_INC="<path_to_netcdf_module>" NETCDF_LIBS="-L<path_to_netcdf_libs> -lnetcdff"
```
or
```
./configure --with-emc EMC_FCFLAGS="-I<path_to_netcdf_module>" EMC_LIBS="-L<path_to_netcdf_libs> -lnetcdff"
```

You need to provide the corresponding path to the netCDF fortran installation module file (`netcdf.mod`).
For example:
```
./configure --with-netcdf NETCDF_INC=/opt/local/include NETCDF_LIBS="-L/opt/local/lib -lnetcdff"
```


## Model setup

Download the EMC model file and put it into this folder `DATA/IRIS_EMC/`.
Create a symbolic link to the downloaded model file with the name `model.nc`, for example like:
```
cd DATA/IRIS_EMC/
ln -s Alaska.JointInversion-RF+Vph+HV-1.Berg.2020-nc4.nc model.nc
```
The mesher will look for a file name `DATA/IRIS_EMC/model.nc` to read in the model values.


## Simulation setup

Then, select in the `DATA/Par_file`:
```
MODEL            =  EMC_model
```
to read in and use the EMC model. 

Additionally, for local EMC models the mesher `xmeshfem3D` can use parameters to produce local meshes, with a cut-off depth and local doubling layers.
In `Par_file`, the following settings could be defined, for example in the Alaska model case:
```
## regional mesh cut-off
# using this flag will cut-off the mesh in the mantle at a layer matching to the given cut-off depth.
# this flag only has an effect for regional simulations, i.e., for NCHUNKS values less than 6.
REGIONAL_MESH_CUTOFF            = .true.

# regional mesh cut-off depth (in km)
# possible selections are: 24.4d0, 80.d0, 220.d0, 400.d0, 600.d0, 670.d0, 771.d0
REGIONAL_MESH_CUTOFF_DEPTH      = 200.d0

# regional mesh cut-off w/ a second doubling layer below 220km interface
# (by default, a first doubling layer will be added below the Moho, and a second one below the 771km-depth layer.
#  Setting this flag to .true., will move the second one below the 220km-depth layer for regional mesh cut-offs only.)
REGIONAL_MESH_ADD_2ND_DOUBLING  = .false.

# flag to turn on local mesh layout
USE_LOCAL_MESH                  = .true.

# total number of mesh layers for local mesh
# (moho used will be the fictitious moho depth, i.e., at 40 or 35 km depth depending on RMOHO_STRETCH_ADJUSTMENT)
NUMBER_OF_LAYERS_CRUST          = 4
NUMBER_OF_LAYERS_MANTLE         = 8

# number of doubling layers
NDOUBLINGS                      = 2
# position of doubling layer (counted from top down)
NZ_DOUBLING_1                   = 2
NZ_DOUBLING_2                   = 5
NZ_DOUBLING_3                   = 0
NZ_DOUBLING_4                   = 0
NZ_DOUBLING_5                   = 0

# time step size
DT                              = 0.05
```
Note that the time step size `DT` can be specified directly in the `Par_file` to overwrite the automatic time step size estimation.
For choosing an appropriate time step, the output in the `output_mesher.txt` file contains information about a "Maximum suggested time step",
which can be used as a guideline to select a stable `DT` size.


## Implementation

Currently, we support the above mentioned Alaska and CVM model files only. 
Other models might define different parameters and are likely not compatible at the moment.

More general support will hopefully come in future.


