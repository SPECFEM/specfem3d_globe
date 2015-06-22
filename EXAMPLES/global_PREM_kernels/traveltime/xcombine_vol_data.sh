#!/bin/bash
#
#
# note: script requires executable 'mesh2vtu'
##############################################
# USER PARAMETERS

# partitions
numCPUs=24

# kernel directory
dir=DATABASES_MPI/

# low (0) / high (1) resolution
res=$1

# colorbar maximum value
maxcolor=2e-07

# combines all partitions
do_all=0

##############################################

if [ "$res" == "" ]; then echo "usage: ./xcombine_vol_data.sh res[0=low/1=high-resolution/2=mid-resolution]"; exit 1; fi

echo
echo "plotting sensitivity kernels (for crust/mantle region)"
echo

# creates slice file
slice=slices_all.txt
if [ "$do_all" == "1" ]; then
  echo "1" | awk '{for(i=0;i<numCPUs;i++)print i}' numCPUs=$numCPUs > $slice
else
  # only for slices on minor arc: source (slice 4), receiver (slice 7)
  echo "4" > $slice
  echo "5" >> $slice
  echo "7" >> $slice
fi

# for visualization
cp -v ~/SPECFEM3D_GLOBE/utils/Visualization/VTK_ParaView/AVS_continent_boundaries.inp .
cp -v ~/SPECFEM3D_GLOBE/utils/Visualization/VTK_ParaView/paraviewpython-example.py .

# only for crust_mantle region
echo
echo "alpha kernel"
echo
./bin/xcombine_vol_data $slice alpha_kernel $dir $dir OUTPUT_FILES/ $res 1
echo ""
mesh2vtu OUTPUT_FILES/reg_1_alpha_kernel.mesh OUTPUT_FILES/reg_1_alpha_kernel.vtu
echo ""
rm -f OUTPUT_FILES/reg_*alpha*.mesh

# colorbar
# see RGBPoints section
echo "maximum color value = $maxcolor"
echo
sed "s:1e-09:$maxcolor:g" state_alpha_kernel.pvsm > tmp_alpha.pvsm

./paraviewpython-example.py tmp_alpha.pvsm
mv -v image.jpg image_alpha_kernel.jpg

echo
echo "beta kernel"
echo
./bin/xcombine_vol_data $slice beta_kernel $dir $dir OUTPUT_FILES/ $res 1
echo ""
mesh2vtu OUTPUT_FILES/reg_1_beta_kernel.mesh OUTPUT_FILES/reg_1_beta_kernel.vtu
echo ""
rm -f OUTPUT_FILES/reg_*beta*.mesh

sed "s:alpha:beta:g" tmp_alpha.pvsm > tmp_beta.pvsm
./paraviewpython-example.py tmp_beta.pvsm
mv -v image.jpg image_beta_kernel.jpg

echo
echo "rho kernel"
echo
./bin/xcombine_vol_data $slice rho_kernel $dir $dir OUTPUT_FILES/ $res 1
echo ""
mesh2vtu OUTPUT_FILES/reg_1_rho_kernel.mesh OUTPUT_FILES/reg_1_rho_kernel.vtu
echo ""
rm -f OUTPUT_FILES/reg_*rho*.mesh

sed "s:alpha:rho:g" tmp_alpha.pvsm > tmp_rho.pvsm
./paraviewpython-example.py tmp_rho.pvsm
mv -v image.jpg image_rho_kernel.jpg

echo
echo "visualize kernel vtu files in directory: OUTPUT_FILES/ using e.g. Paraview"
echo "done"

