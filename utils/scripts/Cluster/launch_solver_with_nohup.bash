#!/bin/bash

# useful for instance on the cluster at LMA Marseille

# put this in file "mymachines" there:
# gpu-1 slots=4 max_slots=8
# gpu-2 slots=4 max_slots=8
# gpu-3 slots=4 max_slots=8

nohup mpirun -np 24 --byslot -machinefile mymachines ./bin/xspecfem3D > output_night.txt 2> output_err.txt &

