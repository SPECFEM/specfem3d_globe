#!/bin/bash

### Dimitri Komatitsch, CNRS Marseille, France, May 2014

# The first plotting command needs and _only_ needs "-K"
# All following but the last need both "-K" and "-O"
# The last one needs _only_ "-O"

if [[ ! -f GMTglobe.cpt ]]; then
    echo "Color palette file GMTglobe.cpt not found!"
    exit 1
fi

ps=topography_of_SPECFEM_mesh.ps

central_meridian=0

rm -f $ps

gmtset ANNOT_FONT_SIZE_PRIMARY 10p HEADER_FONT_SIZE 18p PLOT_DEGREE_FORMAT ddd:mm:ssF

# convert the non-regular cubed-sphere grid output of SPECFEM to a regular grid for GMT
 surface observation_grid_long_lat_topo_for_GMT.txt -Gtopography_of_SPECFEM_mesh.grd -Rd -I4m -f0x,1y -V

## DK DK offset of 4.5cm
## DK DK -Rd is an alias for -R-180/180/-90/90
grdimage topography_of_SPECFEM_mesh.grd -CGMTglobe.cpt -Rd -JK$central_meridian/9i -Y4.5c -K -V > $ps

## DK DK offset of -1.5 inch
pscoast -Rd -JK$central_meridian/9i -B45g30:."Topography of SPECFEM3D_GLOBE mesh": -W -Dc -A1000 -U/-0.75i/-1.5i/"SPECFEM3D_GLOBE calculations by Dimitri Komatitsch" -V -K -O >> $ps

psscale -CGMTglobe.cpt -D12.5/-1.5/16/0.25h -B2000.:"Elevation (m)": -V -O >> $ps

## DK DK convert the final file to PDF
ps2pdf $ps

# Clean up
rm -f .gmt* gmt.conf gmt.history topography_of_SPECFEM_mesh.grd $ps

