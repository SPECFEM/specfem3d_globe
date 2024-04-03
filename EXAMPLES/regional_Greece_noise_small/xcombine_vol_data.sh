#!/bin/bash
#
#
##############################################

# kernel directory
dir="DATABASES_MPI/"

# low (0) / high (1) resolution
res=$1

##############################################

if [ "$res" == "" ]; then echo "usage: ./xcombine_vol_data.sh res[0=low/1=high-resolution/2=mid-resolution]"; exit 1; fi

# for visualization
ln -s ../../utils/Visualization/VTK_ParaView/AVS_continent_boundaries.inp
ln -s ../../utils/Visualization/VTK_ParaView/paraviewpython-example.py

mkdir -p OUTPUT_FILES
if [ ! -f OUTPUT_FILES/sr.vtk ] && [ -f OUTPUT_FILES_1/sr.vtk ]; then cp -v OUTPUT_FILES_1/sr.vtk OUTPUT_FILES/; fi

echo
echo "alpha_kernel"
echo
./bin/xcombine_vol_data_vtk all alpha_kernel $dir $dir OUTPUT_FILES/ $res 1 | tee tmp.log

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

min=`grep "min/max" tmp.log | awk '{print $3 }' | sort | head -n 1`
max=`grep "min/max" tmp.log | awk '{print $4 }' | sort | tail -n 1`
echo
echo "statistics:"
echo "  alpha_kernel min/max: $min $max"
echo

./paraviewpython-example.py state_alpha_kernel.pvsm
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

mv -v image.jpg OUTPUT_FILES/image_alpha_kernel.jpg

echo
echo "beta_kernel"
echo
# only for crust_mantle region
./bin/xcombine_vol_data_vtk all beta_kernel $dir $dir OUTPUT_FILES/ $res 1 | tee tmp.log

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

min=`grep "min/max" tmp.log | awk '{print $3 }' | sort | head -n 1`
max=`grep "min/max" tmp.log | awk '{print $4 }' | sort | tail -n 1`
echo
echo "statistics:"
echo "  beta_kernel min/max: $min $max"
echo

sed "s:alpha:beta:g" state_alpha_kernel.pvsm > tmp_beta.pvsm
./paraviewpython-example.py tmp_beta.pvsm
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

mv -v image.jpg OUTPUT_FILES/image_beta_kernel.jpg

echo
echo "rho_kernel"
echo
./bin/xcombine_vol_data_vtk all rho_kernel $dir $dir OUTPUT_FILES/ $res 1 | tee tmp.log

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

min=`grep "min/max" tmp.log | awk '{print $3 }' | sort | head -n 1`
max=`grep "min/max" tmp.log | awk '{print $4 }' | sort | tail -n 1`
echo
echo "statistics:"
echo "  rho_kernel min/max: $min $max"
echo

sed "s:alpha:rho:g" state_alpha_kernel.pvsm > tmp_rho.pvsm
./paraviewpython-example.py tmp_rho.pvsm
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

mv -v image.jpg OUTPUT_FILES/image_rho_kernel.jpg

echo
echo "visualize kernel vtu files in directory: OUTPUT_FILES/ using e.g. Paraview"
echo "done"

