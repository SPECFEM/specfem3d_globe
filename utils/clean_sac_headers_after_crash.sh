#!/bin/bash

# For this script to work you need to have the sismoutil package available a the
# ORFEUS web-page http://www.orfeus-eu.org

NPTS=`grep 'Total number of time steps written:' OUTPUT_FILES/output_solver.txt | tail -1 | awk '{print $7}'`

echo ' Changing header variable NPTS to ' $NPTS ' ...'

for i in $*
do
  echo " Working on file $i"
  file_ext=`echo $i | sed 's/.*\.\([^\.]*\)/\1/'`
  if [ $file_ext == 'sacan' ]
  then
    vi -e -s -c "16s/ [0-9]*$/\ $NPTS/" -c 'x' $i
  elif [ $file_ext == 'sac' ]
  then
    chsac NPTS $NPTS -f $*
  else
    echo " Cannot clean header variable NPTS of file $i!"
    echo " NOT a SAC file extension!"
  fi
done

echo ' ... done!'

