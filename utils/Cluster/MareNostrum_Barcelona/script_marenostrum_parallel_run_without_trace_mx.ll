#! /bin/ksh

### DK DK submit this with "mnsubmit name_of_script.ll"

#@ job_name = Specfem3D_MPI

#@ initialdir = .

#####
# One node of MareNostrum has two dual-core processors
# therefore the maximum number of tasks per node is four.
#####

#####################################################################
## Running the job with tracing step
#####################################################################
#@ total_tasks = 64
#@ tasks_per_node = 4
############## Wall clock limit hhh:mm:ss
#@ wall_clock_limit = 02:00:00
#@ output = Specfem3D_run_%j.out
#@ error = Specfem3D_run_%j.err
#@ queue
#@ features = mx

#environment
MP_EUILIB=mx
OBJECT_MODE=64
MP_RSH=ssh

    srun ./xgenerate_databases
#    srun ./bin/xspecfem3D

