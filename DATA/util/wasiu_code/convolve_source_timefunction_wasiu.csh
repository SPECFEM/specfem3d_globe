#!/bin/csh

set hdur = 0.04

foreach file ( $* )

set nlines = `wc -l $file `
echo $nlines > input_convolve_code.txt
echo $hdur >> input_convolve_code.txt

echo convolving $file with hdur = $hdur using lines $nlines 

./xconvolve_source_timefunction < $file > ${file}.convolved
rm input_convolve_code.txt

end

