#!/usr/local/bin/gnuplot -persist
#
#    

set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'OUTPUT_FILES/PAS.TS.LHZ.sem.ascii' w l
