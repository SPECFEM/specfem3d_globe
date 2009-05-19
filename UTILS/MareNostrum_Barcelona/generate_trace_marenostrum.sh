#!/bin/sh

export MPITRACE_HOME=/gpfs/apps/CEPBATOOLS/mpitrace-mx/64
export MPTRACE_CONFIG_FILE=mpitrace.xml

export LD_PRELOAD=${MPITRACE_HOME}/lib/libmpitrace.so

## Run the desired program
$*

