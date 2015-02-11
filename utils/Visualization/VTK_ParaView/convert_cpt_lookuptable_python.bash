#!/bin/bash
#
# based on GMT makecpt
#
if [ $# -eq 0 ]; then
  echo "convert_cpt_lookuptable.bash cpt_name_from_gmt (no_green,seis,etc)"
  exit
fi

cpt=$1
ncolors=25
echo "colortable.SetNumberOfTableValues($ncolors)"
makecpt -C$cpt -I -T0/$ncolors/1 > temp.cpt
nline=`wc temp.cpt | awk '{print $1}'`
awk 'NR > 3 && NR < '$nline'-2 {print "colortable.SetTableValue( ",NR-4,",",$2/255,",",$3/255,",",$4/255,",",1.0,")"}' temp.cpt
