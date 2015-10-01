#!/bin/bash
#
# GMT plotting script
#
# plots cross-section output from tool xcreate_cross_section using GMT commands
# (based on GMT version 5.1.2)
#
# usage example:
# - horizontal cross-section
#   > ./plot_cross_section_gmt.sh OUTPUT_FILES/cross_section_vpv.dat 0
#
# - vertical cross-section
#   > ./plot_cross_section_gmt.sh OUTPUT_FILES/cross_section_vpv.dat 1
#
# please modify this script to fit your needs...
#
###############################

# region
R_h=-Rd                   # horizontal
R_v=-R0/360/3500/6371     # vertical

# projection
J_h=-JQ0/8i               # horizontal
J_v=-JP7i                 # vertical

# annotation
B_h=-B45a90/10a45WeSn:."Cross-Section":
B_v=-B15a15/100a100WeSn:."Cross-Section":

# interpolate surface otherwise plot points
use_surface_interpolation=0

# interpolation
I=-I1.0/1.0

# little perspective view
use_perspective=1

################################

# filename
file=$1

# flag for vertical (1) / horizontal (0) cross-section
plot_vert=$2

if [ "$2" == "" ]; then echo "usage: ./plot_cross_section_gmt.sh cross-section-file (e.g. cross_section_vpv.dat) flag-vertical[0==horizontal/1==vertical]"; exit 1; fi

# data
if [ "$plot_vert" == "1" ]; then
  echo ""
  echo "vertical cross-section"
  echo ""
  # extracts data (lon/lat/radius/value/..)
  awk '{if (index($0,"#") == 0){print $1,$2,$3,$4,0.1;}}' $file > $file.gmt
  #awk '{if (index($0,"#") == 0){if ($7 < 1.0) print $1,$2,$3,$4,0.1;}}' $file > $file.gmt

  # data info
  gmtinfo $file.gmt

  lon_min=`gmtinfo -El0 $file.gmt | awk '{print $1}'`
  lon_max=`gmtinfo -Eh0 $file.gmt | awk '{print $1}'`
  lat_min=`gmtinfo -El1 $file.gmt | awk '{print $2}'`
  lat_max=`gmtinfo -Eh1 $file.gmt | awk '{print $2}'`
  r_min=`gmtinfo -El2 $file.gmt | awk '{print $3}'`
  r_max=`gmtinfo -Eh2 $file.gmt | awk '{print $3}'`
  val_min=`gmtinfo -El3 $file.gmt | awk '{print $4}'`
  val_max=`gmtinfo -Eh3 $file.gmt | awk '{print $4}'`

  echo ""
  echo "data statistics:"
  echo "  lon min/max = $lon_min / $lon_max"
  echo "  lat min/max = $lat_min / $lat_max"
  echo "  r   min/max = $r_min   / $r_max"
  echo "  data value min/max = $val_min / $val_max"
  echo
else
  echo ""
  echo "horizontal cross-section"
  echo ""
  # extracts data (lon/lat/value/..)
  awk '{if (index($0,"#") == 0){print $1,$2,$4,0.1;}}' $file > $file.gmt

  # data info
  gmtinfo $file.gmt

  lon_min=`gmtinfo -El0 $file.gmt | awk '{print $1}'`
  lon_max=`gmtinfo -Eh0 $file.gmt | awk '{print $1}'`
  lat_min=`gmtinfo -El1 $file.gmt | awk '{print $2}'`
  lat_max=`gmtinfo -Eh1 $file.gmt | awk '{print $2}'`
  val_min=`gmtinfo -El2 $file.gmt | awk '{print $3}'`
  val_max=`gmtinfo -Eh2 $file.gmt | awk '{print $3}'`

  echo ""
  echo "data statistics:"
  echo "  lon min/max = $lon_min / $lon_max"
  echo "  lat min/max = $lat_min / $lat_max"
  echo "  data value min/max = $val_min / $val_max"
fi

ps_file=$file.ps
grdfile=$file.grd

# project vertical cross section onto plane
if [ "$plot_vert" == "1" ]; then
  echo "vertical projection..."
  echo
  R=$R_v
  J=$J_v
  B=$B_v

  # project vertical cross section onto plane
  # project coordinates on (theta,r) plane
  project $file.gmt -C$lon_min/$lat_min -E$lon_max/$lat_max > $file.gmt2
  awk '{print $6,$3,$4,0.1;}' $file.gmt2 > $file.gmt3

  # determines range for polar plot
  min=`gmtinfo -El0 $file.gmt3 | awk '{print $1}'`
  max=`gmtinfo -Eh0 $file.gmt3 | awk '{print $1}'`
  epi=`echo "scale=0; ($max - $min) / 1" | bc`
  echo "epicentral distance: $epi"

  # adapt plotting range
  if [ $epi -le 355 ]; then
    # adds margin
    theta_min=`echo "scale=0; ($min - 2 ) / 1" | bc`
    theta_max=`echo "scale=0; ($max + 2 ) / 1" | bc`
  else
    theta_min=$min
    theta_max=$max
  fi
  R=-R$theta_min/$theta_max/$r_min/$r_max
  echo "region range: $R"
  echo
else
  echo "horizontal plot..."
  echo
  R=$R_h
  J=$J_h
  B=$B_h
fi

# color scale
#T=-T-0.001/0.001/0.00001
min=`echo "scale=1; ($val_min - 0.5 ) / 1" | bc`
max=`echo "scale=1; ($val_max + 0.5 ) / 1" | bc`
din=`echo "scale=1; ($val_max - $val_min) / 10" | bc`
T=-T$min/$max/$din
echo "color range: $T"
makecpt -Cseis $T -Z -V > $file.cpt

# get start and end points for projection
# first point
lon1=$lon_min
lat1=$lat_min
# second point
lon2=$lon_max
lat2=$lat_max
echo
echo "point locations (lon/lat): ($lon1/$lat1) to ($lon2/$lat2)"
echo


# GMT
# defaults
gmtset MAP_FRAME_TYPE fancy
gmtset PS_MEDIA 10ix10i
gmtset FONT_ANNOT_PRIMARY 9p,Helvetica FONT_TITLE 15p,Helvetica

# starts ps-file
psxy $J $R -K -P -T -V > $ps_file

# global perspective view
if [ "$use_perspective" == "1" ]; then
  echo
  echo "perspective view..."
  echo
  # projection from infinity to lon/lat of first point
  PROJ=-JG$lon1/$lat1/2i
  # positioning
  if [ "$plot_vert" == "1" ]; then
    off_x=2.5i
    off_y=2.5i
  else
    off_x=0i
    off_y=4i
  fi
  pscoast $PROJ -Rg -B5g15 -Glightbrown -Slightblue -W -Dl -P -K -O -V -X$off_x -Y$off_y >> $ps_file

  if [ "$plot_vert" == "1" ]; then
    # vertical cross-section
    # great-circle line
    psxy -J -R -W2 -L -K -O -V <<EOF  >> $ps_file
    $lon1 $lat1
    $lon2 $lat2
EOF
    # circle
    psxy -J -R -Sc0.3 -W1 -G100/100/100 -K -O -V <<EOF  >> $ps_file
    $lon1 $lat1
    $lon2 $lat2
EOF
    # text annotation
    pstext -J -R -K -O -V <<EOF >> $ps_file
    $lon1 $lat1 15 0 1 RT A\030
    $lon2 $lat2 15 0 1 RT B\030
EOF
  else
    # horizontal cross-section
    # region
    psxy -J -R -W2 -L -K -O -V <<EOF  >> $ps_file
    $lon1 $lat1
    $lon1 $lat2
    $lon2 $lat2
    $lon2 $lat1
EOF
  fi
  # shifts back
  psxy -J -R -K -O -X-$off_x -Y-$off_y <<EOF >> $ps_file
EOF
  echo "done perspective"
  echo
fi

# base map
psbasemap $R $J $B -K -O -V -P  >> $ps_file

# plot data points
if [ "$plot_vert" == "1" ]; then
  # simple dot plot
  psxy $file.gmt3 $R $J -C$file.cpt -Sc $B -K -O -V -P >> $ps_file

  # text annotation
  echo "$lon1 $lat1 6371.0 15 0 1 RT A\030" | project -C$lon_min/$lat_min -E$lon_max/$lat_max -Fpz | pstext -J -R -K -O -V >> $ps_file
  echo "$lon2 $lat2 6371.0 15 0 1 RT B\030" | project -C$lon_min/$lat_min -E$lon_max/$lat_max -Fpz | pstext -J -R -K -O -V >> $ps_file

else
  if [ "$use_surface_interpolation" == "1" ]; then
    echo "using surface interpolation..."
    echo
    # grid sampling
    blockmean $R $I $file.gmt -V > $file.gmt.mean
    surface $R $I $file.gmt.mean -G$grdfile -V
    #xyz2grd $file.mean $I $R -G$grdfile -V
    #sphinterpolate $R $I $file.mean -G$grdfile -V
    # info
    grdinfo -L2 $grdfile
    # re-sampling
    grdsample $grdfile -G$grdfile.2 $I -r -V
    grdinfo -L2 $grdfile.2
    # color palette
    # instead of -T0/15/1
    range=`grdinfo $grdfile.2 -T1`
    echo ""
    echo "color range: $range"
    echo ""
    # image plot
    grdimage $grdfile.2 $J $R -C$file.cpt $B -K -V -P > $ps_file
    # contour line
    contour_value=1
    grdcontour $grdfile.2 -J -B -C$contour_value -A5 -Gd3i -S4 -K -O -V >> $ps_file
  else
    # simple dot plot
    psxy $file.gmt $R $J -C$file.cpt -Sc $B -K -O -V -P >> $ps_file
  fi

  # coast line
  pscoast -J -R -W0.1 -Dh -A1000 -K -O -V >> $ps_file
fi

# scale
da=`echo "scale=0; ($val_max - $val_min) / 3" | bc`
echo "scale increment: $da"
echo
#Bscale=-Ba0.001:'':
Bscale=-Ba$da
psscale -D3/-0.5/3/0.2h $Bscale -C$file.cpt -K -O -V >> $ps_file

# ends ps-file
psxy -J -R -O -T -V >> $ps_file

# check
if [[ $? -ne 0 ]]; then echo ""; echo "gmt plotting failed, please check..."; echo ""; exit 1; fi

# convert: ImageMagick command (see http://www.imagemagick.org)
# converting to pdf
convert -quality 100 $ps_file $file.pdf

echo
echo "see file: $file.pdf"
echo

