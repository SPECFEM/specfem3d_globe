#!/bin/sh

###########################################################################
# Shared library preloading example (README)
#
# So as to trace with the dynamic library, just set LD_PRELOAD to the .so
# MPItrace library and set the required environment variables. In this
# example we set MPITRACE_ON to activate the tracing and the MPTRACE_COUNTERS
# to select which counters we wanna obtain from PAPI.
#
###########################################################################

export MPITRACE_HOME=/gpfs/apps/CEPBATOOLS/64.hwc
export LD_PRELOAD=${MPITRACE_HOME}/lib/libmpitrace.so
export MPITRACE_ON=1
# the list of available counters is in /gpfs/apps/CEPBATOOLS/tracing-example/counters.txt
# monitor L1 below
#export MPTRACE_COUNTERS="0x80000000,0x80000017,0x80000032,0x8000003b"
# monitor L2 below
export MPTRACE_COUNTERS="0x80000002,0x80000019,0x80000032,0x8000003b"
# 0x80000000 is PAPI_L1_DCM   # cannot do both L1 and L2 at the same time, need to start two runs
# 0x80000017 is PAPI_L1_LDM
# 0x80000002 is PAPI_L2_DCM
# 0x80000019 is PAPI_L2_LDM
# 0x80000032 is PAPI_TOT_INS
# 0x8000003b is PAPI_TOT_CYC
#
export MPITRACE_MPI_COUNTERS_ON=1
export MPTRACE_COUNTERS_DOMAIN=all
export MPITRACE_MPI_CALLER=1,2,3
#########export MPTRACE_CONFIG_FILE=mpitrace_extended.xml
export MPTRACE_BUFFER_SIZE=150000

# last year 2006: begin
#export MPTRACE_DIR=/gpfs/scratch/hpce07/hpce07084/TRACE_DIR/
#export MPTRACE_FILE_SIZE=10
#export MPTRACE_BUFFER_SIZE=25000
#export MPTRACE_COUNTERS="0x80000004,0x8000000e,0x80000032,0x8000003b"
# last year 2006: end

## Run the desired program
$*

