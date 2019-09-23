#!/bin/bash
#
#
# note: script requires executable 'mesh2vtu'
##############################################
# USER PARAMETERS

# partitions
numCPUs=4

# kernel directory
dir=DATABASES_MPI/

# low (0) / high (1) resolution
res=$1

##############################################

if [ "$res" == "" ]; then echo "usage: ./xcombine_vol_data.sh res[0=low/1=high-resolution/2=mid-resolution]"; exit 1; fi

echo
echo "plotting sensitivity kernels"
echo

# creates slice file
slice=slices_all.txt
echo "1" | awk '{for(i=0;i<numCPUs;i++)print i}' numCPUs=$numCPUs > $slice

# for visualization
cp ../../utils/Visualization/VTK_ParaView/AVS_continent_boundaries.inp ./

kernels=( alpha_kernel \
          beta_kernel \
          rho_kernel )

count=0
for kernel in ${kernels[@]};
do

echo
echo "$kernel"
echo
./bin/xcombine_vol_data $slice $kernel $dir $dir OUTPUT_FILES/ $res 1 > tmp.log
mesh2vtu.pl -i OUTPUT_FILES/reg_1_$kernel.mesh -o OUTPUT_FILES/reg_1_$kernel.vtu >> tmp.log

#mesh2vtu.pl -i OUTPUT_FILES/reg_2_$kernel.mesh -o OUTPUT_FILES/reg_2_$kernel.vtu >> tmp.log
#mesh2vtu.pl -i OUTPUT_FILES/reg_3_$kernel.mesh -o OUTPUT_FILES/reg_3_$kernel.vtu >> tmp.log
rm -f OUTPUT_FILES/reg_*$kernel.mesh

min=`grep "min/max" tmp.log | awk '{print $3 }' | sort | head -n 1`
max=`grep "min/max" tmp.log | awk '{print $4 }' | sort | tail -n 1`
echo "  $kernel min/max: $min $max"

sed "s:alpha_kernel:$kernel:g" state_alpha_kernel.pvsm > tmp.pvsm
./paraviewpython-example.py tmp.pvsm
mv image.jpg image_$kernel.jpg

count=$((count+1))
done

echo
echo "visualize kernel vtu files in directory: OUTPUT_FILES/ using e.g. Paraview"
echo "done"


