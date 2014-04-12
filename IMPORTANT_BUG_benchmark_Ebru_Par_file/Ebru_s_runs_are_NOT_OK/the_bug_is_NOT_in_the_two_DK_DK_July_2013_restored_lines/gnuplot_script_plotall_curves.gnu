#set term x11
set term wxt

set xrange [500:3615]

plot "specfem3d_globe_ifort_July_2013_restored/AAE.IU.MXE.sem.ascii" w l lc 1, "specfem3d_globe_ifort_July_2013_restored_mcmodel_medium/AAE.IU.MXE.sem.ascii" w l lc 3
pause -1 "hit key..."

#############################################

