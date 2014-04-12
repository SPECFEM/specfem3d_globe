#set term x11
set term wxt

set xrange [500:3615]

plot "specfem3d_globe_ifort_mcmodel_medium_s40RTS/AAE.IU.MXE.sem.ascii" w l lc 1, "specfem3d_globe_ifort_s40RTS/AAE.IU.MXE.sem.ascii" w l lc 3, "specfem3d_globe_gfortran4_6_s40RTS/AAE.IU.MXE.sem.ascii" w l lc 4, "specfem3d_globe_portland_pgf90_s40RTS/AAE.IU.MXE.sem.ascii" w l lc 5, "specfem3d_globe_ifort_s40RTS_with_new_Th_Gilet_crustal_hash_table/AAE.IU.MXE.sem.ascii" w l lc 6

pause -1 "hit key..."

