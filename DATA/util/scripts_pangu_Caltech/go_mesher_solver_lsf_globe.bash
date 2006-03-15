#!/bin/bash
#BSUB -o OUTPUT_FILES/%J.o
#BSUB -a mpich_gm
#BSUB -J globe_mesher_solver

if [ -z $USER ]; then
      echo "could not run go_mesher_solver_...bash as no USER env is set"
      exit 2
fi

BASEMPIDIR=/scratch/$USER/DATABASES_MPI

# script to run the mesher

# read DATA/Par_file to get information about the run

# compute total number of nodes needed
NPROC_XI=`grep NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep NCHUNKS DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))

rm -r -f OUTPUT_FILES
mkdir OUTPUT_FILES

# obtain lsf job information
echo "$LSB_MCPU_HOSTS" > OUTPUT_FILES/lsf_machines
echo "$LSB_JOBID" > OUTPUT_FILES/jobid

remap_lsf_machines.pl OUTPUT_FILES/lsf_machines >OUTPUT_FILES/machines

# now just make the dir (cleanup should be done afterwards after collect seismo, otherwise  we clean up another runs seismos)
shmux -M50 -Sall -c "mkdir -p $BASEMPIDIR" - < OUTPUT_FILES/machines >/dev/null

# now modify the Par_file to have the proper path specific to this job id
# this ASSUMES that the Par_file just has a vanilla LOCAL_PATH line
sed -e "s:^LOCAL_PATH .*:LOCAL_PATH                      =  $BASEMPIDIR:" < DATA/Par_file > DATA/Par_file.tmp
mv DATA/Par_file.tmp DATA/Par_file


# not sure why we have to make this on the local node
my_local_path=`grep LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

echo starting MPI mesher on $numnodes processors
echo " "
echo starting run in current directory $PWD
echo " "
echo mesh files will be saved in directory $my_local_path on the nodes
echo " "

#### use this on LSF
mpirun.lsf --gm-no-shmem --gm-copy-env $PWD/xmeshfem3D
mpirun.lsf --gm-no-shmem --gm-copy-env $PWD/xspecfem3D
