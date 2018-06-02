#!/bin/bash
#################################################

# reference adjoint station
network="II"
station="RAYN"
compE="MXE"
compN="MXN"
compZ="MXZ"
en="sem.ascii"

# window start/end time (in seconds)
t_start=30.0
t_end=74.0

#################################################

echo
echo "creating adjoint source for station: $network.$station"
echo

# adjoint sources will be in folder SEM/
currentdir=`pwd`
mkdir -p SEM

# needs traces
sta=$network.$station
if [ ! -e OUTPUT_FILES/$sta.$compE.$en ]; then echo "please make sure traces OUTPUT_FILES/$sta.*.$en are available"; exit 1; fi

rm -f SEM/$sta.*
cp -v OUTPUT_FILES/$sta.*.$en SEM/
echo

# compile adjoint_source tool
if [ ! -e xcreate_adjsrc_traveltime ]; then
  # creates adjoint sources
  cd ../../utils/adjoint_sources/traveltime

  # fortran compiler (as specified in Makefile)
  FC=`grep '^FC .*' ../../../Makefile | cut -d = -f 2 | sed "s/^[ \t]*//"`
  if [ "$FC" == "" ]; then echo "fortran compiler not found, exiting..."; exit 1; fi
  CC=`grep '^CC .*' ../../../Makefile | cut -d = -f 2 | sed "s/^[ \t]*//"`
  if [ "$CC" == "" ]; then echo "C compiler not found, exiting..."; exit 1; fi

  echo "compiling xcreate_adjsrc_traveltime:"
  echo "  using fortran compiler = $FC"
  echo "  using C compiler       = $CC"
  echo

  cp Makefile Makefile.host
  sed -i "s:F90 .*:F90 = $FC:" Makefile.host
  sed -i "s:CC .*:CC = $CC:" Makefile.host

  rm -rf xcreate_adjsrc_traveltime
  make -f Makefile.host
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi

  cp -v xcreate_adjsrc_traveltime $currentdir/SEM/
  cd $currentdir
fi
if [ ! -e SEM/xcreate_adjsrc_traveltime ]; then echo "please make xcreate_adjsrc_traveltime and copy to SEM/"; exit 1; fi

echo
echo "running adjoint source creation"
echo
# creates adjoint sources
cd SEM/

# uses transversal component for adjoint source
type=4
# back-azimuth (station RAYN, measured clockwise from north, from seismic station towards source)
baz=314.6431

./xcreate_adjsrc_traveltime $t_start $t_end $type $sta.*.$en $baz
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

if [ ! -e $sta.$compE.$en.adj ]; then echo "error creating adjoint sources, please check..."; exit 1; fi
echo

mv -v $sta.$compE.$en.adj $sta.$compE.adj
mv -v $sta.$compN.$en.adj $sta.$compN.adj
mv -v $sta.$compZ.$en.adj $sta.$compZ.adj

# create STATIONS_ADJOINT file with adjoint source location
fgrep $station ../DATA/STATIONS > ./STATIONS_ADJOINT
cp -v ./STATIONS_ADJOINT ../DATA/

cd ../

echo
echo

