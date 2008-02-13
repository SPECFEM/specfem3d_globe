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

# @ account = dimitri
# @ class = benchmark

#@ initialdir = .

#####
# One node of MareNostrum has two dual-core processors
# therefore the maximum number of tasks per node is four.
#####

#####################################################################
## Running the job with tracing step
#####################################################################
#@ total_tasks = 2166
#@ tasks_per_node = 4
# Wall clock limit hhh:mm:ss
######### hpce is a queue that has a 24-hour (i.e. 86400 seconds) wallclock time limit
#@ wall_clock_limit = 72:00:00
#@ output = Specfem3D_run_%j.out
#@ error = Specfem3D_run_%j.err
#@ queue

#environment
MP_EUILIB=gm
OBJECT_MODE=64
MP_RSH=ssh

    srun ./xmeshfem3D
    srun ./xspecfem3D

