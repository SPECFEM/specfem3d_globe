#!/bin/bash
#
# The script compares seismograms in two different directories 
#
##################################################
echo
echo "running seismogram comparisons:"
echo

# compares seismograms by plotting correlations
../../utils/compare_seismogram_correlations.py OUTPUT_FILES/ OUTPUT_FILES_reference_OK/

echo
echo "done"
echo
