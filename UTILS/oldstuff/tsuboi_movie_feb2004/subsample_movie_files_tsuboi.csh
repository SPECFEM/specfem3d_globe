#!/bin/csh

# Intel Linux f90 compiler
ifort -fast -o xsubsample_movie_files_tsuboi subsample_movie_files_tsuboi.f90

foreach file ( $* )

echo processing file $file

./xsubsample_movie_files_tsuboi < $file > new_$file

end

rm -r xsubsample_movie_files_tsuboi

