#!/bin/bash
######################################
# USER PARAMETERS

# number of time steps
nsteps=$1

# time step size
dt=$2

# period range
Tmin=2
Tmax=100

#####################################

if [ "$1" == "" -o "$2" == "" ]; then echo "usage: ./xgenerate_noise_source.sh nsteps dt"; exit 1; fi

# creates noise spectrum (using new low noise model)
matlab -nodisplay << EOF
NOISE_TOMOGRAPHY($nsteps,$dt,$Tmin,$Tmax,'NLNM')
exit
EOF
cp -v S_squared NOISE_TOMOGRAPHY/


