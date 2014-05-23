#!/bin/bash

# The first plotting command needs - and _only_ needs "-K"
# All following but the last need both, "-K" and "-O"
# The last one needs _only_ "-O"

central_meridian=0

gmtset ANNOT_FONT_SIZE_PRIMARY 10p HEADER_FONT_SIZE 18p PLOT_DEGREE_FORMAT ddd:mm:ssF

# DK DK for absolute norm of the g vector
#makecpt  -T9.0/9.2/0.001 -Z > color2.cpt
# DK DK for the difference between the 3D calculation and the 1D reference calculation
#makecpt  -T-0.03/0.03/0.0001 -Z > color2.cpt
#-0.00668716430664 z_max: 0.00217151641846 name: z
makecpt  -T-0.007/0.007/0.00001 -Z > color2.cpt

################################

file=results_norm_of_g

## DK DK subtract the 3D result with topography and s20rts from the 1D result for the reference ellipsoid calculation
grdmath new_SPECFEM3D_GLOBE_NEX_160_on_600_cores_whole_Earth_with_topo_and_s40rts/${file}.grd new_SPECFEM3D_GLOBE_NEX_160_on_600_cores_whole_Earth_with_ellipticity_only_and_PREM1D/${file}.grd NEG ADD = ${file}_difference.grd

## DK DK offset de 4.5cm
## DK DK -Rd signifie -R-180/180/-90/90
grdimage ${file}_difference.grd  -Ccolor2.cpt -Rd -JK$central_meridian/9i -Y4.5c -K -V > ${file}_difference.ps

## DK DK offset de -1.5 inch
pscoast -Rd -JK$central_meridian/9i -B45g30:."norm_of_g_difference": -W -Dc -A1000 -U/-0.75i/-1.5i/"SPECFEM3D_GLOBE gravity calculations by Dimitri Komatitsch" -V -O -K >> ${file}_difference.ps

psscale -Ccolor2.cpt -D12.5/-1.5/16/0.25h -B0.006:"Difference with reference ellipsoid (in m/s2)": -V -O >> ${file}_difference.ps

## DK DK convert the final file to PDF
ps2pdf ${file}_difference.ps
rm -f ${file}_difference.ps

# Clean up
rm -f .gmt* gmt.conf gmt.history color2.cpt

