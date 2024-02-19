#!/bin/bash

\rm -f output_night.o output_night.e core.*

ulimit -S -s unlimited

ccc_msub -p gen6351 script_MPI_150_solver_CEA_CCRT.sh

