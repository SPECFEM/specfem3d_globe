#!/bin/bash

# this is optional, the default topography file topo_bathy_etopo4_smoothed_window_7.dat used by SPECFEM3D_GLOBE
# is already included in the source code

# this only needs to be done once (and for all)

cd DATA/topo_bathy

echo " "
echo "downloading the whole topography database, which has a size of 837 MB; this may take a while..."
echo " "

wget http://data.geodynamics.org/specfem/specfem3d_globe/topo_bathy/etopo1_ice_c.grd.bz2
wget http://data.geodynamics.org/specfem/specfem3d_globe/topo_bathy/etopo1_ice_c_resampled_at_2minutes.grd.bz2
wget http://data.geodynamics.org/specfem/specfem3d_globe/topo_bathy/etopo2v2c.grd.bz2
wget http://data.geodynamics.org/specfem/specfem3d_globe/topo_bathy/topo_bathy_etopo1_ice_c_original_unmodified_unsmoothed.bin.bz2
wget http://data.geodynamics.org/specfem/specfem3d_globe/topo_bathy/topo_bathy_etopo1_ice_c_resampled_at_2minutes_original_unmodified_unsmoothed.bin.bz2
wget http://data.geodynamics.org/specfem/specfem3d_globe/topo_bathy/topo_bathy_etopo1_ice_c_resampled_at_2minutes_smoothed_window_3.bin.bz2
wget http://data.geodynamics.org/specfem/specfem3d_globe/topo_bathy/topo_bathy_etopo1_ice_c_smoothed_window_3.bin.bz2
wget http://data.geodynamics.org/specfem/specfem3d_globe/topo_bathy/topo_bathy_etopo2v2c_original_unmodified_unsmoothed.bin.bz2
wget http://data.geodynamics.org/specfem/specfem3d_globe/topo_bathy/topo_bathy_etopo2v2c_smoothed_window_3.bin.bz2
wget http://data.geodynamics.org/specfem/specfem3d_globe/topo_bathy/topo_bathy_etopo4_from_etopo2_subsampled.bin.bz2

echo " "
echo "You can now uncompress the topography file you want using the bunzip2 command."
echo " "
