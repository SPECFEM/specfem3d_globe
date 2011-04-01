#!/bin/bash -eu

DIR_LOCAL=$1
DIR_GLOBAL=$2
iproc=`printf "%06d" $PBS_VNODENUM`

##### other files can be collected in the same way
mv $DIR_LOCAL/*$iproc*reg1*kernel*.bin      $DIR_GLOBAL/
mv $DIR_LOCAL/*$iproc*reg1_*array_dims.txt  $DIR_GLOBAL/
mv $DIR_LOCAL/*$iproc*reg1_solver_data*.bin $DIR_GLOBAL/

