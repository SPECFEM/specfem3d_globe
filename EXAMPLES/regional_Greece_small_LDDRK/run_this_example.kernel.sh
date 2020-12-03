#!/bin/bash
#
# regional simulation example to compute an adjoint kernel
#
# script runs mesher and solver for an adjoint kernel using mpirun
# on 4 CPUs
#
# modify accordingly for your own system specifics
##################################################

echo "running example: `date`"
currentdir=`pwd`

echo
echo "  setting up example..."
echo
./setup_this_example.kernel.bash
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# run mesher & solver
echo
echo "  running script..."
echo
./run_mesher_solver.kernel.bash
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo `date`

