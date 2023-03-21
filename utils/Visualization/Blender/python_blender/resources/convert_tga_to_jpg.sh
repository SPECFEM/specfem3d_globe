#!/bin/bash


if [ "$1" == "" ]; then echo "usage: ./convert_tga_to_jpg.sh <filename>"; exit 1; fi

file=$1

## tga
#convert resources/cloud_combined_8192.tif -compress none -alpha off -colorspace sRGB -channel rgb -type truecolor clouds.tga
#identify clouds.tga


## ppm
# not working, tif is rgb with 8-bit 0-255 values
# tif -> ppm
#file=resources/cloud_combined_8192.tif
#gdalinfo -mm $file
# directly output ppm file in 16-bit
#gdal_translate -of PNM -ot uint16 -scale -10898 8271 0 65535 -r bilinear -outsize 8192 0 $file topo_8192.ppm
# should be grayscale 16-bit
#identify -verbose topo_8192.ppm


# jpg ending

filename=`basename $file`

ext="${filename##*.}"
name="${filename%.*}"

echo "file: $file"
echo
echo "name     : $name"
echo "extension: $ext"
echo

identify $file
identify -format "Colors: %k\nFile Size: %b\n" $file
echo

# imagemagick: only converts to 8-bit images
# not working...
# 16-bit gray to 16-bit jpeg
#options="-auto-level -depth 12" # -quality 100%
options="-quality 100%"

if [[ "$name" == *topo* ]]; then
# topo ppm file has 16-bit, uses png format to have again 16-bit grayscale
out=${name}.png
options="-type grayscale"
else
out=${name}.jpg
fi



echo "output file: $out"
echo

# convert to jpg with imageMagick
convert $options $file $out
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo

# might need flipping
if [ "$name" == "world" ]||[ "$name" == "night" ]||[ "$name" == "earth" ]||[ "$name" == "clouds" ]; then
echo "flipping image..."
convert $out -flip $out
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi


fi



identify $out
identify -format "Colors: %k\nFile Size: %b\n" $out
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo


