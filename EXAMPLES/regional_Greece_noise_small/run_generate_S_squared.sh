#!/bin/bash
######################################
# USER PARAMETERS
#
# e.g. call by: ./run_generate_S_squared.sh 3599  0.169376865
#
# number of time steps
noise_nsteps=$1

# time step size
dt=$2

# period range (simulation minimum resolution ~15s)
Tmin=10  #2
Tmax=30  #100

# noise model (new low noise model)
noise_model=NLNM

#####################################

if [ "$1" == "" -o "$2" == "" ]; then echo "usage: ./run_generate_S_squared.sh nsteps dt [Tmin] [Tmax]"; exit 1; fi

# optional
if [ "$3" != "" ]; then Tmin=$3; fi
if [ "$4" != "" ]; then Tmax=$4; fi

echo "simulation setup:"
echo "  DT    = $dt"
echo "  Tmin/Tmax = $Tmin / $Tmax"
echo
echo "noise:"
echo "  number of steps = $noise_nsteps"
echo "  model           = $noise_model"
echo

mkdir -p NOISE_TOMOGRAPHY/

# creates noise spectrum (using new low noise model)
# generates S_squared file
if [ 0 == 1 ]; then
  ## matlab
  matlab -nosplash -nodisplay -r "NOISE_TOMOGRAPHY($noise_nsteps,$dt,$Tmin,$Tmax,\'$noise_model\');exit"
else
  ## python & numpy
  ./NOISE_TOMOGRAPHY.py $noise_nsteps $dt $Tmin $Tmax $noise_model
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# move to noise setup folder
mv -v S_squared NOISE_TOMOGRAPHY/

echo
echo "done"
echo

