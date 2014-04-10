
#set term x11
set term wxt

set xrange [500:3615]

plot "specfem3d_globe_ifort_no_ftz/AAE.IU.MXE.sem.ascii" w l lc 1, "specfem3d_globe_ifort/AAE.IU.MXE.sem.ascii" w l lc 7, "specfem3d_globe_ifort_mcmodel_medium_for_solver_only/AAE.IU.MXE.sem.ascii" w l lc 11
pause -1 "hit key..."

plot "specfem3d_globe_ifort_no_mcmodel_medium_with_mcmodel_medium_created_mesh/AAE.IU.MXE.sem.ascii" w l lc 2, "specfem3d_globe_GNU_4_8/AAE.IU.MXE.sem.ascii" w l lc 3, "specfem3d_globe_mcmodel_medium_GNU_4_8/AAE.IU.MXE.sem.ascii" w l lc 5, "specfem3d_globe_ifort_mcmodel_medium_no_ftz/AAE.IU.MXE.sem.ascii" w l lc 6, "specfem3d_globe_ifort_mcmodel_medium/AAE.IU.MXE.sem.ascii" w l lc 9, "specfem3d_globe_ifort_mcmodel_medium_for_mesher_only/AAE.IU.MXE.sem.ascii" w l lc 10
pause -1 "hit key..."

plot "specfem3d_globe_mcmodel_medium_GNU_4_6/AAE.IU.MXE.sem.ascii" w l lc 4, "specfem3d_globe_GNU_4_6/AAE.IU.MXE.sem.ascii" w l lc 8
pause -1 "hit key..."

#############################################

plot "specfem3d_globe_ifort_no_ftz/AAE.IU.MXE.sem.ascii" w l lc 1, "specfem3d_globe_ifort_no_mcmodel_medium_with_mcmodel_medium_created_mesh/AAE.IU.MXE.sem.ascii" w l lc 2, "specfem3d_globe_mcmodel_medium_GNU_4_6/AAE.IU.MXE.sem.ascii" w l lc 4
pause -1 "hit key..."

#############################################

#plot "specfem3d_globe_ifort_no_ftz/AAE.IU.MXE.sem.ascii" w l lc 1, "specfem3d_globe_ifort_no_mcmodel_medium_with_mcmodel_medium_created_mesh/AAE.IU.MXE.sem.ascii" w l lc 2, "specfem3d_globe_GNU_4_8/AAE.IU.MXE.sem.ascii" w l lc 3, "specfem3d_globe_mcmodel_medium_GNU_4_6/AAE.IU.MXE.sem.ascii" w l lc 4, "specfem3d_globe_mcmodel_medium_GNU_4_8/AAE.IU.MXE.sem.ascii" w l lc 5, "specfem3d_globe_ifort_mcmodel_medium_no_ftz/AAE.IU.MXE.sem.ascii" w l lc 6, "specfem3d_globe_ifort/AAE.IU.MXE.sem.ascii" w l lc 7, "specfem3d_globe_GNU_4_6/AAE.IU.MXE.sem.ascii" w l lc 8, "specfem3d_globe_ifort_mcmodel_medium/AAE.IU.MXE.sem.ascii" w l lc 9
#pause -1 "hit key..."

