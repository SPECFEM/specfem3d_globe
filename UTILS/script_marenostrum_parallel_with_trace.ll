#! /bin/ksh

#@ job_name = Specfem3D_MPI

### type "llclass" to see the standard queues
### NAME            MAX_CPUS           wall_clock_time
###################################################
### interactive         1                     2h
### debug              64                    10m
### hpce              512                    48h
###################################################
######### hpce is a slow queue with 512 processors
#@ class = hpce

#@ group = hpce07

#@ initialdir = .

#####
# One node of MareNostrum has two dual-core processors
# therefore the maximum number of tasks per node is four.
#####

#####################################################################
## Running the job with tracing step
#####################################################################
#@ step_name = trace_step
#@ total_tasks = 96
# Wall clock limit hhh:mm:ss
######### hpce is a queue that has a 24-hour (i.e. 86400 seconds) wallclock time limit
#@ wall_clock_limit = 02:00:00
#@ job_type = parallel
#@ blocking  = unlimited
#@ node_usage = not_shared
#@ restart = no
#@ output = Specfem3D_run_$(jobid).out
#@ error = Specfem3D_run_$(jobid).err
#@ queue

#####################################################################
## Parallel merging step (using 3 tasks)
#####################################################################
#@ dependency = (trace_step == 0)
#@ step_name = merge_step
#@ output = Specfem3D_merge_$(jobid).out
#@ error = Specfem3D_merge_$(jobid).err
#@ total_tasks = 3
#@ wall_clock_limit = 01:00:00    # Modify this if needed
#@ job_type = parallel
#@ blocking  = unlimited
#@ node_usage = not_shared
#@ restart = no
#@ queue
#####################################################################

MPITRACE_HOME=/gpfs/apps/CEPBATOOLS/64.hwc

#environment
MP_EUILIB=gm
OBJECT_MODE=64
MP_RSH=ssh

# we get the machine list and number of processors from LoadLeveler
MLIST=machine_list.$$.${LOADL_STEP_NAME}
/opt/ibmll/LoadL/full/bin/ll_get_machine_list > ${MLIST}
NPROCS=`cat ${MLIST} | wc -l`

case ${LOADL_STEP_NAME} in

      trace_step)
    mpirun -np ${NPROCS} -machinefile $LL_MACHINE_LIST ./xmeshfem3D
    mpirun -np ${NPROCS} -machinefile $LL_MACHINE_LIST ./generate_trace_marenostrum.sh ./xspecfem3D
    ;;

      merge_step)
    mpirun -np ${NPROCS} -machinefile ${MLIST} ${MPITRACE_HOME}/bin/mpimpi2prv -f TRACE.mpits -maxmem 1024 -syn -o trace.prv
    ;;

      *)
    echo "Uknown step ${LOADL_STEP_NAME}"
    ;;

esac

