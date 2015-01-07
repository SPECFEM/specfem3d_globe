#set term postscript eps monochrome dashed "Helvetica" 22

set term x11
 
set size 1,0.65
 
set xlabel "Time (s)"
set ylabel "Displacement (cm)"
 
set nozeroaxis

set xrange [400:3000]
set yrange [-0.002:0.002]

#set output "bolivia_PAS_longit.eps"
plot "PAS_longit_semd_final.txt" t 'SEM longitudinal from 2002 GJI paper' w l lc 1, "PAS_longit_modes_final.txt" t 'Modes' w l lc 3, "CI.PAS.MXR.sem.ascii.convolved" t 'SEM results given by your run' w l lc 6
pause -1 "hit key"

#set output "bolivia_PAS_transv.eps"
plot "PAS_transv_semd_final.txt" t 'SEM transverse from 2002 GJI paper' w l lc 1, "PAS_transv_modes_final.txt" t 'Modes' w l lc 3, "CI.PAS.MXT.sem.ascii.convolved" t 'SEM results given by your run' w l lc 6
pause -1 "hit key"

#set output "bolivia_PAS_vertical.eps"
plot "PAS_vert_semd_final.txt" t 'SEM vertical from 2002 GJI paper' w l lc 1, "PAS_vert_modes_final.txt" t 'Modes' w l lc 3, "CI.PAS.MXZ.sem.ascii.convolved" t 'SEM results given by your run' w l lc 6
pause -1 "hit key"

