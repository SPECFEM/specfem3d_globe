#!/bin/csh

# add logo banner to all images before creating a movie

foreach file ($*)

echo adding logo to image $file...

set file1 = `basename $file 0.tiff`

set file2 = `echo $file1 | sed 's/^imagemovie0//'`

montage -tile 1x2 -frame 0 -geometry +0+0 -depth 8 -strip $file ../banner_logos.tiff ../sequence_finale/te{$file2}.tif

end

