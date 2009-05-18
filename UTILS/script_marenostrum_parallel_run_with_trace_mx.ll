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
#@ total_tasks = 32
#@ tasks_per_node = 4
# Wall clock limit hhh:mm:ss
######### hpce is a queue that has a 24-hour (i.e. 86400 seconds) wallclock time limit
#@ wall_clock_limit = 00:20:00
#@ output = Specfem3D_run_%j.out
#@ error = Specfem3D_run_%j.err
#@ queue
#@ features = mx

# DK DK this below OK in May 2009
MPITRACE_HOME=/gpfs/apps/CEPBATOOLS/mpitrace-mx/64

#environment
MP_EUILIB=mx
OBJECT_MODE=64
MP_RSH=ssh

    srun ./generate_trace_paraver.sh ./xspecfem2D

# then merge the trace at the end
    sleep 5
    mnsubmit script_marenostrum_parallel_merge_trace.ll

