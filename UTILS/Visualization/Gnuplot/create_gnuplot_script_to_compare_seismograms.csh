#!/bin/csh

#echo set term postscript color solid "Helvetica" 22
#echo set output \"all_seismograms_comparison.ps\"

echo set term pdf color solid
echo set output \"all_seismograms_comparison.pdf\"

#echo "set term x11"
#echo "set term wxt"

echo set xrange \[0:3300\]

foreach file ( OUTPUT_FILES/*.sem.ascii )

  set newfile = `basename $file`

  echo plot \"OUTPUT_FILES_reference_OK/$newfile\" w l lc 1, \"$file\" w l lc 3
# use the line below instead if you have an old install of Gnuplot (the Gnuplot syntax was different back then)
# echo plot \"OUTPUT_FILES_reference_OK/$newfile\" w l 1, \"$file\" w l 3

# uncomment this only when outputting to the screen (X11 or wxt)
#  echo "pause -1 'hit any key...'"

end

