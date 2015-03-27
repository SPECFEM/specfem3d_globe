#!/bin/bash
#PBS -N NOISE_sim
#PBS -q tromp
#PBS -o zzz_global_long/job_info/output
#PBS -e zzz_global_long/job_info/error
#PBS -l nodes=75:ppn=8,walltime=100:00:00

DIR_RESULT="zzz_global_long"
NSTEP=31599

cd $PBS_O_WORKDIR
######################################## pre-simulation ###########################################
## gather simulation infomation
NPROC_XI=`grep ^NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep ^NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep ^NCHUNKS DATA/Par_file | cut -d = -f 2 `
numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))

cat  $PBS_NODEFILE > $DIR_RESULT/job_info/nodes
echo $PBS_JOBID    > $DIR_RESULT/job_info/id

## check number of procs
NPROC_1=`wc -l < $DIR_RESULT/job_info/nodes`  # number of procs assigned to the job
NPROC_2=$numnodes                             # number of procs needed in SPECFEM
if [ $NPROC_1 -ne $NPROC_2 ]; then
   echo "inconsistent number of procs, check your job submission script & DATA/Par_file" 1>&2
   echo "NPROC from job submission: $NPROC_1" 1>&2
   echo "NPROC from DATA/Par_file: $NPROC_2" 1>&2
   exit
fi

## clean local nodes
LOCAL_PATH=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 | sed 's/ //g'`

d=`date`
echo "Start cleaning local disk on nodes, $d"
pbsdsh $PBS_O_WORKDIR/run_NOISE_clean_create.sh   $LOCAL_PATH
d=`date`
echo "Finish cleaning local disk on nodes, $d"

## clean OUTPUT_FILES
rm -rf   OUTPUT_FILES
mkdir -p OUTPUT_FILES

######################################## simulation ###############################################
################################################################### step 1: ensemble forward source
cp DATA/Par_file_NOISE_1_attenuation DATA/Par_file
mpiexec $PBS_O_WORKDIR/bin/xmeshfem3D_attenuation
mpiexec $PBS_O_WORKDIR/bin/xspecfem3D_attenuation

echo "2" > temp_input_movie; echo "1" >> temp_input_movie; echo "$NSTEP" >> temp_input_movie
$PBS_O_WORKDIR/bin/xcreate_movie_AVS_DX_attenuation < temp_input_movie

# backup step 1 results (with attenuation)
DIR_TEMP=$DIR_RESULT/step_1
rm   -rf $DIR_TEMP
mkdir -p $DIR_TEMP
mv OUTPUT_FILES/*   $DIR_TEMP/
cp -v DATA/Par_file    $DIR_TEMP/
cp -v DATA/CMTSOLUTION $DIR_TEMP/
cp -v DATA/STATIONS    $DIR_TEMP/

################################################################ step 2: ensemble forward wavefield
cp DATA/Par_file_NOISE_2_attenuation DATA/Par_file
mpiexec $PBS_O_WORKDIR/bin/xmeshfem3D_attenuation
mpiexec $PBS_O_WORKDIR/bin/xspecfem3D_attenuation

echo "2" > temp_input_movie; echo "1" >> temp_input_movie; echo "$NSTEP" >> temp_input_movie
$PBS_O_WORKDIR/xcreate_movie_AVS_DX_attenuation < temp_input_movie

# backup step 2 results (with attenuation)
DIR_TEMP=$DIR_RESULT/step_2
rm   -rf $DIR_TEMP
mkdir -p $DIR_TEMP
mv OUTPUT_FILES/*   $DIR_TEMP/
cp -v DATA/Par_file    $DIR_TEMP/
cp -v DATA/CMTSOLUTION $DIR_TEMP/
cp -v DATA/STATIONS    $DIR_TEMP/

############################################################## prepare adjoint sources ############
cp -v $DIR_RESULT/step_2/A7.II.MXZ.sem.ascii ./SEM/
./NOISE_adj

################################################################### step 1: ensemble forward source
cp DATA/Par_file_NOISE_1_noattenuation DATA/Par_file
mpiexec $PBS_O_WORKDIR/xmeshfem3D_noattenuation
mpiexec $PBS_O_WORKDIR/xspecfem3D_noattenuation

echo "2" > temp_input_movie; echo "1" >> temp_input_movie; echo "$NSTEP" >> temp_input_movie
$PBS_O_WORKDIR/xcreate_movie_AVS_DX_noattenuation < temp_input_movie

# backup step 1 results (without attenuation)
DIR_TEMP=$DIR_RESULT/step_1_no_attenuation
rm   -rf $DIR_TEMP
mkdir -p $DIR_TEMP
mv OUTPUT_FILES/*   $DIR_TEMP/
cp -v DATA/Par_file    $DIR_TEMP/
cp -v DATA/CMTSOLUTION $DIR_TEMP/
cp -v DATA/STATIONS    $DIR_TEMP/

################################################################ step 2: ensemble forward wavefield
cp DATA/Par_file_NOISE_2_noattenuation DATA/Par_file
mpiexec $PBS_O_WORKDIR/xmeshfem3D_noattenuation
mpiexec $PBS_O_WORKDIR/xspecfem3D_noattenuation

echo "2" > temp_input_movie; echo "1" >> temp_input_movie; echo "$NSTEP" >> temp_input_movie
$PBS_O_WORKDIR/xcreate_movie_AVS_DX_noattenuation < temp_input_movie

# backup step 2 results (without attenuation)
DIR_TEMP=$DIR_RESULT/step_2_no_attenuation
rm   -rf $DIR_TEMP
mkdir -p $DIR_TEMP
mv OUTPUT_FILES/*   $DIR_TEMP/
cp -v DATA/Par_file    $DIR_TEMP/
cp -v DATA/CMTSOLUTION $DIR_TEMP/
cp -v DATA/STATIONS    $DIR_TEMP/

############################################################## step 3: ensemble kernels calculation
cp DATA/Par_file_NOISE_3_noattenuation DATA/Par_file
mpiexec $PBS_O_WORKDIR/xmeshfem3D_noattenuation
mpiexec $PBS_O_WORKDIR/xspecfem3D_noattenuation

echo "2" > temp_input_movie; echo "1" >> temp_input_movie; echo "$NSTEP" >> temp_input_movie
$PBS_O_WORKDIR/xcreate_movie_AVS_DX_noattenuation < temp_input_movie

## collect outputs (e.g., kernels, meshes) from local nodes to global disk
d=`date`
echo "Start collecting outputs from nodes, $d"
pbsdsh $PBS_O_WORKDIR/run_NOISE_collect.sh  $LOCAL_PATH      $PBS_O_WORKDIR/OUTPUT_FILES/
d=`date`
echo "Finish collecting outputs from nodes, $d"

## combine_vol_data
echo "0" > temp_input_movie
for (( iproc=1;iproc<=$numnodes-1;iproc++ ))
do
  echo "$iproc" >> temp_input_movie
done
$PBS_O_WORKDIR/xcombine_vol_data temp_input_movie alpha_kernel $PBS_O_WORKDIR/OUTPUT_FILES $PBS_O_WORKDIR/OUTPUT_FILES $PBS_O_WORKDIR/OUTPUT_FILES 0 1
rm temp_input_movie

## backup step 3 results (without attenuation)
DIR_TEMP=$DIR_RESULT/step_3_no_attenuation
rm   -rf $DIR_TEMP
mkdir -p $DIR_TEMP
mv OUTPUT_FILES/*   $DIR_TEMP/
cp -v DATA/Par_file    $DIR_TEMP/
cp -v DATA/CMTSOLUTION $DIR_TEMP/
cp -v DATA/STATIONS    $DIR_TEMP/

## clean local nodes
pbsdsh $PBS_O_WORKDIR/run_NOISE_clean_create.sh   $LOCAL_PATH
d=`date`
echo "Finish cleaning local disk on nodes, $d"

