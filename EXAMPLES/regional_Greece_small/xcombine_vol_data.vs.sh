#!/bin/bash
#
#
# note: script requires executable 'mesh2vtu'
##############################################
# partitions
numCPUs=4
# slice file
echo "1" | awk '{for(i=0;i<numCPUs;i++)print i}' numCPUs=$numCPUs > slices_all.txt
slice="slices_all.txt"
# kernel directory
dir="DATABASES_MPI/"
# low (0) / high (1) resolution
res="0"

# for visualization
cp ../../UTILS/Visualization/Paraview/AVS_continent_boundaries.inp ./

# velocity model for: vs
par="vs"
echo
echo "velocity model: $par"
echo
./bin/xcombine_vol_data $slice $par $dir $dir OUTPUT_FILES/ $res > tmp.log
mesh2vtu.pl -i OUTPUT_FILES/reg_1_$par.mesh -o OUTPUT_FILES/reg_1_$par.vtu >> tmp.log
mesh2vtu.pl -i OUTPUT_FILES/reg_2_$par.mesh -o OUTPUT_FILES/reg_2_$par.vtu >> tmp.log
mesh2vtu.pl -i OUTPUT_FILES/reg_3_$par.mesh -o OUTPUT_FILES/reg_3_$par.vtu >> tmp.log
rm -f OUTPUT_FILES/reg_*$par*.mesh
min=`grep "min/max" tmp.log | awk '{print $3 }' | sort | head -n 1`
max=`grep "min/max" tmp.log | awk '{print $4 }' | sort | tail -n 1`
echo "  $par min/max: $min $max"
./paraviewpython-example.py state_vs.pvsm
mv image.jpg image_vs.jpg

echo
echo "visualize mesh vtu files in directory: OUTPUT_FILES/ using e.g. Paraview"
echo "done"


