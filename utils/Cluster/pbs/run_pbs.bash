#!/bin/bash
# use the normal queue unless otherwise directed

#module load openmpi/intel-11.1

rm -f OUTPUT_FILES/*

d=`date`
echo "Starting compilation $d"
make clean
make xmeshfem3D
make xcreate_header_file
./bin/xcreate_header_file
make xspecfem3D
d=`date`
echo "Finished compilation $d"

echo "Submitting job"

# mesher
qsub < go_mesher_sge.bash

# solver
#qsub < go_solver_sge.bash
