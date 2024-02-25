#!/bin/bash

# If you want other topographic models, you may download the whole topography database

# this is optional, the default topography file topo_bathy_etopo4_smoothed_window_7.dat used by SPECFEM3D_GLOBE
# is already included in the source code

# you will need to call this script from the root directory of the package.
# this only needs to be done once (and for all)

if [ ! -e DATA/topo_bathy ]; then
  echo "you will need to call this script from the root directory of the package..."
  exit 1
fi


cd DATA/topo_bathy

echo " "
echo "downloading the whole topography database, which has a size of 837 MB; this may take a while..."
echo " "

# note: the geodynamics data server is currently offline.
#       topo/bathy data can also be created directly by the python script 'run_create_topo_bathy_file.py' in the DATA/topo_bathy/ directory.
#
#       as alternative, you can also download topo/bathy data from the following github repository:
#       https://github.com/SPECFEM/specfem-data

# deprecated
#url="http://data.geodynamics.org/specfem/specfem3d_globe/topo_bathy/"
#wget ${url}/etopo1_ice_c.grd.bz2
#wget ${url}/etopo1_ice_c_resampled_at_2minutes.grd.bz2
#wget ${url}/etopo2v2c.grd.bz2
#wget ${url}/topo_bathy_etopo1_ice_c_original_unmodified_unsmoothed.bin.bz2
#wget ${url}/topo_bathy_etopo1_ice_c_resampled_at_2minutes_original_unmodified_unsmoothed.bin.bz2
#wget ${url}/topo_bathy_etopo1_ice_c_resampled_at_2minutes_smoothed_window_3.bin.bz2
#wget ${url}/topo_bathy_etopo1_ice_c_smoothed_window_3.bin.bz2
#wget ${url}/topo_bathy_etopo2v2c_original_unmodified_unsmoothed.bin.bz2
#wget ${url}/topo_bathy_etopo2v2c_smoothed_window_3.bin.bz2
#wget ${url}/topo_bathy_etopo4_from_etopo2_subsampled.bin.bz2

# via github.com/SPECFEM/specfem-data
url="https://github.com/SPECFEM/specfem-data/blob/master/topo_bathy/"

# uncomment desired file(s)
files=( \
        ##topo_bathy_etopo1_smoothed_window_3.dat.bin.bz2 \
        ##topo_bathy_etopo1_unmodified_unsmoothed.dat.bin.bz2 \
        topo_bathy_etopo2_smoothed_window_3.dat.bin.bz2 \
        ##topo_bathy_etopo2_unmodified_unsmoothed.dat.bin.bz2 \
        ##topo_bathy_etopo4_smoothed_window_7.dat.bin.bz2 \
        ##topo_bathy_etopo4_unmodified_unsmoothed.dat.bin.bz2 \
       )
for file in ${files[@]};
do
  echo "download file: $file"
  sleep 2
  wget -O $file ${url}/${file}?raw=true
  echo ""
done
echo "downloaded files into directory: DATA/topo_bathy/"

echo " "
echo "You can now uncompress the topography file you want using the bunzip2 command."
echo "You will also need to uncomment it in file setup/constants.h.in before configuring and compiling the code."
echo " "


