#!/bin/bash
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
if [ ! -e ./paraviewpython-example.py ]; then
ln -s ../../utils/Visualization/VTK_ParaView/paraviewpython-example.py
ln -s ../../utils/Visualization/VTK_ParaView/AVS_continent_boundaries.inp
fi

kernels=( alpha_kernel \
          beta_kernel \
          rho_kernel )

count=0
for kernel in ${kernels[@]};
do

echo
echo "kernel: $kernel"
echo
./bin/xcombine_vol_data_vtk $slice $kernel $dir $dir OUTPUT_FILES/ $res 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo

sed "s:alpha_kernel:$kernel:g" state_alpha_kernel.pvsm > tmp.pvsm

./paraviewpython-example.py tmp.pvsm $count


mv -v image$count.jpg OUTPUT_FILES/image$count.$kernel.jpg
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

count=$((count+1))
done

echo
echo "visualize kernel vtu files in directory: OUTPUT_FILES/ using e.g. Paraview"
echo "done"


