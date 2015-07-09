#!/bin/bash

# remove all the old mesh files in the local /scratch of each blade
# on "pangu" at Caltech

# Dimitri Komatitsch, University of Pau, November 2007

if [ -z $USER ]; then
  echo "cannot run this script because no USER env is set"
  exit 2
fi

BASEMPIDIR=/scratch/$USER

echo cleaning local scratch space $BASEMPIDIR on each node of the cluster

grep compute- /opt/lsfhpc/conf/lsf.cluster.lsfhpc | expand | cut -f 1 -d ' ' > ___________bubu

shmux -M 50 -S all -c "rm -r -f $BASEMPIDIR" - < ___________bubu >/dev/null

rm ___________bubu

