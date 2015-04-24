#!/bin/bash
######################################
# USER PARAMETERS
#
# e.g. call by: ./xgenerate_noise_source.sh 3599  0.169376865
#
# number of time steps
nsteps=$1

# time step size
dt=$2

# period range (simulation minimum resolution ~15s)
Tmin=10  #2
Tmax=100 #100

#####################################

if [ "$1" == "" -o "$2" == "" ]; then echo "usage: ./xgenerate_noise_source.sh nsteps dt"; exit 1; fi

# creates noise spectrum (using new low noise model)
matlab -nodisplay << EOF
NOISE_TOMOGRAPHY($nsteps,$dt,$Tmin,$Tmax,'NLNM')
exit
EOF

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# moves to NOISE_TOMOGRAPHY setup directory
echo
mkdir -p NOISE_TOMOGRAPHY
mv -v S_squared NOISE_TOMOGRAPHY/

echo
echo "done"
echo

