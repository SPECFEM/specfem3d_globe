#!/bin/bash
#
##############################################
# USER PARAMETERS

# kernel directory
dir=DATABASES_MPI/

# low (0) / high (1) resolution
res=$1

##############################################

if [ "$res" == "" ]; then echo "usage: ./xcombine_vol_data.sh res[0=low/1=high-resolution/2=mid-resolution]"; exit 1; fi

echo
echo "plotting sensitivity kernels"
echo

# for visualization
ln -s ../../utils/Visualization/VTK_ParaView/AVS_continent_boundaries.inp

kernels=( alpha_kernel \
          beta_kernel \
          rho_kernel )

for kernel in ${kernels[@]};
do

echo
echo "$kernel"
echo
# crust/mantle == region 1
./bin/xcombine_vol_data_vtu all $kernel $dir $dir OUTPUT_FILES/ $res 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

sed "s:alpha_kernel:$kernel:g" state_alpha_kernel.pvsm > tmp.pvsm
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

./paraviewpython-example.py tmp.pvsm
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

mv -v image.jpg image_$kernel.jpg

done

echo
echo "visualize kernel vtu files in directory: OUTPUT_FILES/ using e.g. Paraview"
echo "done"


