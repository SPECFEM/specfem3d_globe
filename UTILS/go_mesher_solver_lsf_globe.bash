#!/bin/bash -v
#BSUB -o OUTPUT_FILES/%J.o
#BSUB -a mpich_gm
#BSUB -J go_mesher_solver_lsf

BASEMPIDIR=/scratch/$USER/DATABASES_MPI

echo "$LSB_MCPU_HOSTS" > OUTPUT_FILES/lsf_machines
echo "$LSB_JOBID" > OUTPUT_FILES/jobid

remap_lsf_machines.pl OUTPUT_FILES/lsf_machines >OUTPUT_FILES/machines

shmux -M50 -Sall -c "rm -r -f /scratch/$USER; mkdir -p /scratch/$USER; mkdir -p $BASEMPIDIR" - < OUTPUT_FILES/machines >/dev/null

current_pwd=$PWD

mpirun.lsf --gm-no-shmem --gm-copy-env $current_pwd/xmeshfem3D

mpirun.lsf --gm-no-shmem --gm-copy-env $current_pwd/xspecfem3D

# collect seismograms and clean up
mkdir -p SEM
cd SEM
collect_seismo.pl ../OUTPUT_FILES/machines
cleanbase.pl ../OUTPUT_FILES/machines

