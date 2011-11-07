#!/bin/csh

# script by Dimitri Komatitsch, CNRS Marseille, November 2011
# to compare two sets of seismograms (for instance to make sure that a modification made in the code does not change the seismograms)

# call it with:
#   ./create_gnuplot_script_to_compare_seismograms.csh OUTPUT_FILES/*.semd > plotall.gnu
# (assuming that the seismogram files are called *.semd in directory OUTPUT_FILES)
# and place the reference seismograms to compare to in directory REFERENCE/OUTPUT_FILES

# also make sure that parameter "hdur" is greater than zero in file DATA/CMTSOLUTION
# because only seismograms convolved with a source time function should be compared
# (seismograms computed for a Heaviside source time function before convolution cannot be reliably compared
# because they are too sensitive to roundoff noise)

#echo set term postscript color solid "Helvetica" 22
#echo set output \"tutu.ps\"

echo set term X11

#echo set xrange \[0:14000\]

foreach file ($*)

# change "w l 1" to "w l lc 1" to use more recent Gnuplot syntax
# (i.e., draw a line between data points, with line color 1)
# (same thing with "w l 3")
echo plot \"REFERENCE/$file\" w l 1, \"$file\" w l 3
echo pause -1 \"Hit any key to see the next seismogram...\"

end

