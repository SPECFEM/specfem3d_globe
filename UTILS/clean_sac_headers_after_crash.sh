#!/bin/bash

# For this script to work you need to have the sismoutil package available a the
# ORFEUS web-page http://www.orfeus-eu.org

NPTS=`grep 'Total number of time steps written' OUTPUT_FILES/output_solver.txt | tail -1 | sed -e 's/.*[:space:]//'`

echo ' Changing header variable NPTS to ' $NPTS

chsac NPTS $NPTS -f $*

