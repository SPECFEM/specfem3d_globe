#!/bin/bash
#
#
# note: script requires executable 'mesh2vtu'
##############################################

# partitions
numCPUs=4

# kernel directory
dir="DATABASES_MPI/"

# low (0) / high (1) resolution
res="1"

##############################################

# slice file
echo "1" | awk '{for(i=0;i<numCPUs;i++)print i}' numCPUs=$numCPUs > slices_all.txt
slice="slices_all.txt"

# for visualization
cp -v ~/SPECFEM3D_GLOBE/utils/Visualization/VTK_ParaView/AVS_continent_boundaries.inp .
cp -v ~/SPECFEM3D_GLOBE/utils/Visualization/VTK_ParaView/paraviewpython-example.py .

echo
echo "alpha_kernel"
echo
./bin/xcombine_vol_data $slice alpha_kernel $dir $dir OUTPUT_FILES/ $res 1 > tmp.log
mesh2vtu.pl -i OUTPUT_FILES/reg_1_alpha_kernel.mesh -o OUTPUT_FILES/reg_1_alpha_kernel.vtu >> tmp.log
rm -f OUTPUT_FILES/reg_*alpha*.mesh
min=`grep "min/max" tmp.log | awk '{print $3 }' | sort | head -n 1`
max=`grep "min/max" tmp.log | awk '{print $4 }' | sort | tail -n 1`
echo "  alpha_kernel min/max: $min $max"

./paraviewpython-example.py state_alpha_kernel.pvsm
mv image.jpg image_alpha_kernel.jpg

echo
echo "beta_kernel"
echo
# only for crust_mantle region
./bin/xcombine_vol_data $slice beta_kernel $dir $dir OUTPUT_FILES/ $res 1 > tmp.log
mesh2vtu.pl -i OUTPUT_FILES/reg_1_beta_kernel.mesh -o OUTPUT_FILES/reg_1_beta_kernel.vtu >> tmp.log
rm -f OUTPUT_FILES/reg_*beta*.mesh
min=`grep "min/max" tmp.log | awk '{print $3 }' | sort | head -n 1`
max=`grep "min/max" tmp.log | awk '{print $4 }' | sort | tail -n 1`
echo "  beta_kernel min/max: $min $max"

sed "s:alpha:beta:g" state_alpha_kernel.pvsm > tmp_beta.pvsm
./paraviewpython-example.py tmp_beta.pvsm
mv image.jpg image_beta_kernel.jpg

echo
echo "rho_kernel"
echo
./bin/xcombine_vol_data $slice rho_kernel $dir $dir OUTPUT_FILES/ $res 1 > tmp.log
mesh2vtu.pl -i OUTPUT_FILES/reg_1_rho_kernel.mesh -o OUTPUT_FILES/reg_1_rho_kernel.vtu >> tmp.log
rm -f OUTPUT_FILES/reg_*rho*.mesh
min=`grep "min/max" tmp.log | awk '{print $3 }' | sort | head -n 1`
max=`grep "min/max" tmp.log | awk '{print $4 }' | sort | tail -n 1`
echo "  rho_kernel min/max: $min $max"

sed "s:alpha:rho:g" state_alpha_kernel.pvsm > tmp_rho.pvsm
./paraviewpython-example.py tmp_rho.pvsm
mv image.jpg image_rho_kernel.jpg

echo
echo "visualize kernel vtu files in directory: OUTPUT_FILES/ using e.g. Paraview"
echo "done"


