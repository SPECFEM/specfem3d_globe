#!/bin/bash

/bin/rm -f output_night*.o output_night*.e core.* seismo*txt

/bin/rm -f -r DATABASES_MPI
/bin/rm -f -r OUTPUT_FILES

mkdir -p DATABASES_MPI OUTPUT_FILES
mkdir -p OUTPUT_FILES/DATABASES_MPI

if [[ -f saved_observation_grid_real_x_y_z_used_by_the_code.txt ]]; then
  cp saved_observation_grid_real_x_y_z_used_by_the_code.txt OUTPUT_FILES
fi

ulimit -S -s unlimited

#ccc_msub -p genb002 -q gpu script_MPI_030.sh
ccc_msub -q standard script_MPI_128.sh

