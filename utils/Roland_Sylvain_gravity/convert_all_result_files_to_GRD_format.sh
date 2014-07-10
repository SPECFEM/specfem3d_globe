#!/bin/bash

# from the GMT manual:
# As a rule of thumb, we use blockmean for most smooth data (such as potential fields)
# and blockmedian (or modes) for rough, non-Gaussian data (such as topography).

 blockmean -Rd -I4m -V results_g_x_for_GMT.txt > ________temp________.xyz
 surface ________temp________.xyz -Gresults_g_x.grd -Rd -I4m -f0x,1y -V

 blockmean -Rd -I4m -V results_g_y_for_GMT.txt > ________temp________.xyz
 surface ________temp________.xyz -Gresults_g_y.grd -Rd -I4m -f0x,1y -V

 blockmean -Rd -I4m -V results_g_z_for_GMT.txt > ________temp________.xyz
 surface ________temp________.xyz -Gresults_g_z.grd -Rd -I4m -f0x,1y -V

 blockmean -Rd -I4m -V results_norm_of_g_for_GMT.txt > ________temp________.xyz
 surface ________temp________.xyz -Gresults_norm_of_g.grd -Rd -I4m -f0x,1y -V

 blockmean -Rd -I4m -V results_G_xx_for_GMT.txt > ________temp________.xyz
 surface ________temp________.xyz -Gresults_G_xx.grd -Rd -I4m -f0x,1y -V

 blockmean -Rd -I4m -V results_G_yy_for_GMT.txt > ________temp________.xyz
 surface ________temp________.xyz -Gresults_G_yy.grd -Rd -I4m -f0x,1y -V

 blockmean -Rd -I4m -V results_G_zz_for_GMT.txt > ________temp________.xyz
 surface ________temp________.xyz -Gresults_G_zz.grd -Rd -I4m -f0x,1y -V

 blockmean -Rd -I4m -V results_G_xy_for_GMT.txt > ________temp________.xyz
 surface ________temp________.xyz -Gresults_G_xy.grd -Rd -I4m -f0x,1y -V

 blockmean -Rd -I4m -V results_G_xz_for_GMT.txt > ________temp________.xyz
 surface ________temp________.xyz -Gresults_G_xz.grd -Rd -I4m -f0x,1y -V

 blockmean -Rd -I4m -V results_G_yz_for_GMT.txt > ________temp________.xyz
 surface ________temp________.xyz -Gresults_G_yz.grd -Rd -I4m -f0x,1y -V

# Clean up
rm -f .gmt* ________temp________.xyz

