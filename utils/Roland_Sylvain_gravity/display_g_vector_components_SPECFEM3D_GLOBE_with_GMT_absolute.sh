#!/bin/bash

# The first plotting command needs - and _only_ needs "-K"
# All following but the last need both, "-K" and "-O"
# The last one needs _only_ "-O"

central_meridian=0

gmtset ANNOT_FONT_SIZE_PRIMARY 10p HEADER_FONT_SIZE 18p PLOT_DEGREE_FORMAT ddd:mm:ssF

# DK DK for absolute norm of the g vector
makecpt -Chaxby -T-10/10/0.001 > color2.cpt

################################

file=results_g_x

## DK DK offset de 4.5cm
## DK DK -Rd signifie -R-180/180/-90/90
grdimage new_SPECFEM3D_GLOBE_NEX_160_on_600_cores_whole_Earth_with_topo_and_s40rts/${file}.grd -Ccolor2.cpt -Rd -JK$central_meridian/9i -Y4.5c -K -V > ${file}_absolute.ps

## DK DK offset de -1.5 inch
pscoast -Rd -JK$central_meridian/9i -B45g30:."g_x_absolute": -W -Dc -A1000 -U/-0.75i/-1.5i/"SPECFEM3D_GLOBE gravity calculations by Dimitri Komatitsch" -V -O -K >> ${file}_absolute.ps

psscale -Ccolor2.cpt -D12.5/-1.5/16/0.25h -B0.006:"Absolute value (in m/s2)": -V -O >> ${file}_absolute.ps

## DK DK convert the final file to PDF
ps2pdf ${file}_absolute.ps
rm -f ${file}_absolute.ps

################################

file=results_g_y

## DK DK offset de 4.5cm
## DK DK -Rd signifie -R-180/180/-90/90
grdimage new_SPECFEM3D_GLOBE_NEX_160_on_600_cores_whole_Earth_with_topo_and_s40rts/${file}.grd -Ccolor2.cpt -Rd -JK$central_meridian/9i -Y4.5c -K -V > ${file}_absolute.ps

## DK DK offset de -1.5 inch
pscoast -Rd -JK$central_meridian/9i -B45g30:."g_y_absolute": -W -Dc -A1000 -U/-0.75i/-1.5i/"SPECFEM3D_GLOBE gravity calculations by Dimitri Komatitsch" -V -O -K >> ${file}_absolute.ps

psscale -Ccolor2.cpt -D12.5/-1.5/16/0.25h -B0.006:"Absolute value (in m/s2)": -V -O >> ${file}_absolute.ps

## DK DK convert the final file to PDF
ps2pdf ${file}_absolute.ps
rm -f ${file}_absolute.ps

################################

file=results_g_z

## DK DK offset de 4.5cm
## DK DK -Rd signifie -R-180/180/-90/90
grdimage new_SPECFEM3D_GLOBE_NEX_160_on_600_cores_whole_Earth_with_topo_and_s40rts/${file}.grd -Ccolor2.cpt -Rd -JK$central_meridian/9i -Y4.5c -K -V > ${file}_absolute.ps

## DK DK offset de -1.5 inch
pscoast -Rd -JK$central_meridian/9i -B45g30:."g_z_absolute": -W -Dc -A1000 -U/-0.75i/-1.5i/"SPECFEM3D_GLOBE gravity calculations by Dimitri Komatitsch" -V -O -K >> ${file}_absolute.ps

psscale -Ccolor2.cpt -D12.5/-1.5/16/0.25h -B0.006:"Absolute value (in m/s2)": -V -O >> ${file}_absolute.ps

## DK DK convert the final file to PDF
ps2pdf ${file}_absolute.ps
rm -f ${file}_absolute.ps

################################

# Clean up
rm -f .gmt* gmt.conf gmt.history color2.cpt

