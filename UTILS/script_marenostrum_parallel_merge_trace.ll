#! /bin/ksh

### DK DK submit this with "mnsubmit name_of_script.ll"

#@ job_name = Specfem3D_MPI

#@ initialdir = .

#####
# One node of MareNostrum has two dual-core processors
# therefore the maximum number of tasks per node is four.
#####

#####################################################################
## parallel merging step
#####################################################################
#@ output = Specfem3D_merge_%j.out
#@ error = Specfem3D_merge_%j.err
#@ total_tasks = 32
#@ tasks_per_node = 4
#@ wall_clock_limit = 00:20:00
#@ queue
#@ features = mx 
#####################################################################

#################MPITRACE_HOME=/gpfs/apps/CEPBATOOLS/64.hwc
MPITRACE_HOME=/gpfs/apps/CEPBATOOLS/mpitrace/64
###############MPITRACE_HOME=/gpfs/apps/CEPBATOOLS/64.burst

#environment
MP_EUILIB=mx
OBJECT_MODE=64
MP_RSH=ssh

  srun ${MPITRACE_HOME}/bin/mpimpi2prv -f TRACE.mpits -maxmem 1024 -syn -o trace.prv

