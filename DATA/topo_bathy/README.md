# Topography and bathymetry data

This folder contains binary global topography/bathymetry files for **Earth**, **Mars** and **Moon**. Due to file sizes for higher resolutions, we limit the provided resolutions in this folder to ETOPO4.

**The default resolution is ETOPO4, as defined and set in file `setup/constants.h`. <br>
How can I use higher resolution topography files?**

Higher resolution files can either be created by
- running the script `run_create_topo_bathy_file.py` in this folder

  For example, for creating a binary file with ETOPO1-resolution for **Earth**, type:
  ```
  ./run_create_topo_bathy_file.py etopo1
  ```
  or, for **Mars** with a 2 arc-minute resolution, type:
  ```
  ./run_create_topo_bathy_file.py mars2
  ```
  or, for **Moon** with a 1 arc-minute resolution, type:
  ```
  ./run_create_topo_bathy_file.py moon1
  ```



- downloading pre-created binary files from the data repository [SPECFEM specfem-data](https://github.com/SPECFEM/specfem-data)
  Please check out the data repository for further infos.

You will then need to modify accordingly the corresponding entries in file `setup/constants.h`:

**Earth**
```
!--- Default
! Topography defaults to ETOPO4
  integer, parameter :: EARTH_NX_BATHY = NX_BATHY_4
  integer, parameter :: EARTH_NY_BATHY = NY_BATHY_4
  double precision, parameter :: EARTH_RESOLUTION_TOPO_FILE = RESOLUTION_TOPO_FILE_4
  character (len=*), parameter :: EARTH_PATHNAME_TOPO_FILE = PATHNAME_TOPO_FILE_4
```

**Mars**
```
!--- Default
! Topography defaults to 4-minute
  integer, parameter :: MARS_NX_BATHY = MARS_NX_BATHY_4
  integer, parameter :: MARS_NY_BATHY = MARS_NY_BATHY_4
  double precision, parameter :: MARS_RESOLUTION_TOPO_FILE = MARS_RESOLUTION_TOPO_FILE_4
  character (len=*), parameter :: MARS_PATHNAME_TOPO_FILE = MARS_PATHNAME_TOPO_FILE_4
```

**Moon**
```
!--- Default
! Topography defaults to 4-minute
  integer, parameter :: MOON_NX_BATHY = MOON_NX_BATHY_4
  integer, parameter :: MOON_NY_BATHY = MOON_NY_BATHY_4
  double precision, parameter :: MOON_RESOLUTION_TOPO_FILE = MOON_RESOLUTION_TOPO_FILE_4
  character (len=*), parameter :: MOON_PATHNAME_TOPO_FILE = MOON_PATHNAME_TOPO_FILE_4
```


Note:
To avoid meshing issues (negative Jacobians) with distorted mesh elements due to a rough topography, we usually apply a small smoothing filter to the raw downloaded topo data files. This smoothing step stabilizes the numerical simulation and is chosen empirically. You can test with/without this smoothing filter by modifying the script `run_create_topo_bathy_file.py`.
