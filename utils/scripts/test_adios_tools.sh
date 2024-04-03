#!/bin/bash
#
# tests ADIOS file handling tools
#########################################################

# partitions
NPROC=$1

# paraview plotting
vis=$2

# low (0) / high (1) resolution
#res="1"

#########################################################

if [ "$NPROC" == "" ]; then echo "usage: ./test_adios_tools.sh NPROC[e.g., 4] [plot]"; exit 1; fi


echo
echo "testing tools:"
echo "  NPROC = $NPROC"
if [ "$vis" == "plot" ]; then
  echo "  plotting on"
else
  echo "  plotting off"
fi
echo
echo

# slice file
echo "1" | awk '{for(i=0;i<numCPUs;i++)print i}' numCPUs=$NPROC > slices_all.txt
slice="slices_all.txt"

# for visualization
if [ "$vis" == "plot" ]; then
  cp -v ../../utils/Visualization/VTK_ParaView/AVS_continent_boundaries.inp ./
  cp -v ../../utils/Visualization/VTK_ParaView/paraviewpython-example.py ./
fi

## for kernel summation xsum_kernel
# folder topo/ holds mesh topology (mesh_parameters.bin)
if [ ! -e topo ]; then
ln -s ./DATABASES_MPI/ topo
fi

# folder OUTPUT_SUM/ holds summation result
mkdir -p OUTPUT_SUM
rm -rf OUTPUT_SUM/kernels_sum.bp

# for input kernels
mkdir -p INPUT_KERNELS
# kernel file
if [ ! -e kernels_list.txt ]; then
  echo "event_1" > kernels_list.txt
  cd INPUT_KERNELS/
  if [ ! -e event_1 ]; then
    ln -s ../OUTPUT_FILES/ event_1
  fi
  cd ../
fi

# combines & plots alpha_kl
echo
echo "*****************************************"
echo "xcombine_vol_data_vtk_adios and plot"
echo "*****************************************"
echo
#./bin/xcombine_vol_data_vtk_adios slices_all.txt alpha_kl OUTPUT_FILES/kernels.bp DATABASES_MPI/solver_data.bp OUTPUT_FILES/ 1 1
# checks exit code
#if [[ $? -ne 0 ]]; then exit 1; fi

./bin/xcombine_vol_data_vtk_adios slices_all.txt bulk_betav_kl OUTPUT_FILES/kernels.bp DATABASES_MPI/solver_data.bp OUTPUT_FILES/ 1 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# plotting
if [ "$vis" == "plot" ]; then
  echo
  sed "s:alpha_kl:bulk_betav_kl:g" state_alpha_kl.pvsm > tmp.pvsm
  ./paraviewpython-example.py tmp.pvsm
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
  mv -v image.jpg OUTPUT_FILES/image_betav_kernel-1.jpg
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
fi

mv -v OUTPUT_FILES/reg_1_bulk_betav_kl.vtk OUTPUT_FILES/reg_1_bulk_betav_kl-1.vtk

# sums tiso kernels
echo
echo "*****************************************"
echo "xsum_kernels_adios"
echo "*****************************************"
echo
mpirun -np $NPROC ./bin/xsum_kernels_adios
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# combines bulk_betav_kl from summed kernels
echo
echo "*****************************************"
echo "xcombine_vol_data_vtk_adios w/ summed kernel"
echo "*****************************************"
echo
./bin/xcombine_vol_data_vtk_adios slices_all.txt bulk_betav_kl OUTPUT_SUM/kernels_sum.bp DATABASES_MPI/solver_data.bp OUTPUT_FILES/ 1 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# plotting
if [ "$vis" == "plot" ]; then
  echo
  ./paraviewpython-example.py tmp.pvsm
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
  mv -v image.jpg OUTPUT_FILES/image_betav_kernel-2.jpg
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
fi


# interpolation (onto itself)
echo
echo "*****************************************"
echo "xinterpolate_model_adios"
echo "*****************************************"
echo
mpirun -np $NPROC ./bin/xinterpolate_model_adios DATABASES_MPI/ DATABASES_MPI/ DATABASES_MPI/ OUTPUT_FILES/
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "vsh version: adios"
echo
./bin/xcombine_vol_data_vtk_adios slices_all.txt vsh DATABASES_MPI/model_gll.bp DATABASES_MPI/solver_data.bp OUTPUT_FILES/ 1 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi


echo
echo "vsh version: interpolated"
echo
./bin/xcombine_vol_data_vtk_adios slices_all.txt vsh OUTPUT_FILES/model_gll_interpolated.bp DATABASES_MPI/solver_data.bp OUTPUT_FILES/ 1 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi


# converts to binary files
echo
echo "*****************************************"
echo "xconvert_model_file_adios"
echo "*****************************************"
echo
mpirun -np $NPROC ./bin/xconvert_model_file_adios 1 DATABASES_MPI/ OUTPUT_FILES/
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi


# combines bulk_betav_kl from converted kernels
echo
echo "*****************************************"
echo "xcombine_vol_data_vtk_adios w/ converted model"
echo "*****************************************"
echo
cp -v DATA/Par_file DATA/Par_file.tmp
sed -i "s:^ADIOS_ENABLED .*:ADIOS_ENABLED = .false.:" DATA/Par_file

mpirun -np $NPROC ./bin/xmeshfem3D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# restores original Par_file
mv -v DATA/Par_file.tmp DATA/Par_file

echo
echo "vsh version: adios"
echo
./bin/xcombine_vol_data_vtk_adios slices_all.txt vsh DATABASES_MPI/model_gll.bp DATABASES_MPI/solver_data.bp OUTPUT_FILES/ 1 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "vsh version: binary"
echo
./bin/xcombine_vol_data_vtk slices_all.txt vsh DATABASES_MPI/ OUTPUT_FILES/ OUTPUT_FILES/ 1 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "eta version: adios"
echo
./bin/xcombine_vol_data_vtk_adios slices_all.txt eta DATABASES_MPI/model_gll.bp DATABASES_MPI/solver_data.bp OUTPUT_FILES/ 1 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "eta version: binary"
echo
./bin/xcombine_vol_data_vtk slices_all.txt eta DATABASES_MPI/ OUTPUT_FILES/ OUTPUT_FILES/ 1 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# checks if smoothing is working
echo
echo "*****************************************"
echo "xsmooth_sem_adios"
echo "*****************************************"
echo

GPU_MODE=`grep ^GPU_MODE DATA/Par_file | cut -d = -f 2 `
if [[ "$GPU_MODE" =~ ".true." ]]; then
mpirun -np $NPROC ./bin/xsmooth_sem_adios 1.0 1.0 bulk_betav_kl,bulk_betah_kl OUTPUT_FILES/kernels.bp DATABASES_MPI/solver_data.bp OUTPUT_FILES/kernels_smooth.bp .true.
else
mpirun -np $NPROC ./bin/xsmooth_sem_adios 1.0 1.0 bulk_betav_kl,bulk_betah_kl OUTPUT_FILES/kernels.bp DATABASES_MPI/solver_data.bp OUTPUT_FILES/kernels_smooth.bp
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi


# checks if laplace smoothing is working
echo
echo "*****************************************"
echo "xsmooth_laplacian_sem_adios"
echo "*****************************************"
echo

mpirun -np $NPROC ./bin/xsmooth_laplacian_sem_adios 1.0 1.0 bulk_betav_kl,bulk_betah_kl OUTPUT_FILES/kernels.bp DATABASES_MPI/ OUTPUT_FILES/kernels_smooth_laplace.bp

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo "test successful!"
echo

