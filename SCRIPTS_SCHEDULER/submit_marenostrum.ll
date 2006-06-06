#!/bin/bash

####### LoadLeveler script for MareNostrum in Barcelona (IBM PowerPC970)

#@ job_name = Specfem3D


#### debug is a fast queue with 64 processors and a 10-minute wallclock time limit
######@ class = debug
######@ group = hpce07

#### papi is a fast queue with 64 processors and a 12-hour wallclock time limit
######@ class = papi
######@ group = papi

#### bench is a large queue with 1024 processors but wallclock time must be small
######@ class = bench
######@ group = bench

#### hpce is a slow queue with 512 processors and a 24-hour wallclock time limit
#@ class = hpce
#@ group = hpce07


# Wall clock limit for this job, in hh:mm:ss
#@ wall_clock_limit = 14:00:00

#@ node = 147
#@ tasks_per_node = 2
#@ node_usage = not_shared

#@ initialdir = /home/hpce07/hpce07084/SPECFEM3D_GLOBE

#@ job_type = parallel

#@ error = Specfem3D_$(jobid).err
#@ output = Specfem3D_$(jobid).out

#@ queue

#environment
export MP_EUILIB=gm
export OBJECT_MODE=64
export MP_RSH=ssh

cd /home/hpce07/hpce07084/SPECFEM3D_GLOBE

rm -f machine_list_*
export MLIST=machine_list_$LOADL_STEP_ID
/opt/ibmll/LoadL/full/bin/ll_get_machine_list > $MLIST

NPROCS=`cat $MLIST |wc -l`

rm -r -f OUTPUT_FILES
mkdir OUTPUT_FILES

# create a file containing the jobID
echo $LOADL_STEP_ID | cut -d '.' -f 3 > OUTPUT_FILES/jobid

mpirun -np ${NPROCS} -machinefile $MLIST ./xmeshfem3D

# with trace
#######mpirun -np ${NPROCS} -machinefile $MLIST ./trace ./xspecfem3D

# without trace
mpirun -np ${NPROCS} -machinefile $MLIST ./xspecfem3D

####cd /gpfs/scratch/hpce07/hpce07084/TRACE_DIR/
#####llsubmit /gpfs/scratch/hpce07/hpce07084/TRACE_DIR/trace.cmd

