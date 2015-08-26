#!/bin/bash
#
# The script compares seismograms in two different directories
#
##################################################
echo
echo "running seismogram comparisons:"
echo

# uncompress seismograms
if [ -e OUTPUT_FILES_reference_OK/II.AAK.MXE.sem.ascii.bz2 ]; then
  echo
  echo "unzipping references..."
  echo
  bunzip2 OUTPUT_FILES_reference_OK/*.bz2
  echo
  echo
fi

# compares seismograms by plotting correlations
../../utils/compare_seismogram_correlations.py OUTPUT_FILES/ OUTPUT_FILES_reference_OK/

echo
echo "done"
echo
