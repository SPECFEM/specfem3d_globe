
set term x11

plot "prem.dat" us 1:2 t 'vp' w l 1, "prem.dat" us 1:3 t 'vs' w l 2

#########plot "prem.dat" us 1:2 t 'vp' w l 1, "prem.dat" us 1:3 t 'vs' w l 2, "prem.dat" us 1:4 t 'rho' w l 3

pause -1 "hit key"

