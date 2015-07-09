#!/bin/bash -eu

DIR_LOCAL=$1
iproc=$PBS_VNODENUM

if [ $(($iproc % 8)) = 0 ]; then
   rm -rf   $DIR_LOCAL
   mkdir -p $DIR_LOCAL
fi

