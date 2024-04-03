#!/bin/bash
#
#
##############################################

# kernel directory
dir="DATABASES_MPI/"

# low (0) / high (1) resolution
res=$1

##############################################

if [ "$res" == "" ]; then echo "usage: ./xcombine_vol_data.vsv.sh res[0=low/1=high-resolution/2=mid-resolution]"; exit 1; fi

# for visualization
ln -s ../../utils/Visualization/VTK_ParaView/AVS_continent_boundaries.inp
ln -s ../../utils/Visualization/VTK_ParaView/paraviewpython-example.py

# velocity model for: vs
parameters=( vsv )

for par in ${parameters[@]};
do

echo
echo "velocity model: $par"
echo
./bin/xcombine_vol_data_vtu all $par $dir $dir OUTPUT_FILES/ $res
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

./paraviewpython-example.py state_vsv.pvsm
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

mv -v image.jpg image_$par.jpg

done

echo
echo "visualize mesh vtu files in directory: OUTPUT_FILES/ using e.g. Paraview"
echo "done"


