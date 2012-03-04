#!/bin/csh

# create the four MPEG movie sequences from all the OpenDX *.tif files with logos added separately

#set quality = 100
set quality = 85

 convert -quality $quality te*.tif movie_transparent_europe.mpeg

 convert -quality $quality tp*.tif movie_transparent_pacific.mpeg

 convert -quality $quality oe*.tif movie_opaque_europe.mpeg

 convert -quality $quality op*.tif movie_opaque_pacific.mpeg

