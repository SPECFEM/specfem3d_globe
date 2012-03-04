function [xs,ys,zs,nx,ny,nz] = tp2norm(lats,lons,latr,lonr)
% find out the coordinates of the source (xs,ys,zs) 
% the normal of the plane (nx,ny,nz) 
% constructed by (ths,phs), (thr,phr), and the origin
% for a "banana cross-section"

ts = (90 - lats) * pi/180;
ps = lons * pi / 180;
tr = (90 - latr) * pi/180;
pr = lonr * pi / 180;

[xs,ys,zs] = tp2xyz(ts,ps);
[xr,yr,zr] = tp2xyz(tr,pr);
nx = ys*zr - zs*yr;
ny = -(xs*zr - zs*xr);
nz = xs*yr - ys*xr;
nr = sqrt(nx*nx+ny*ny+nz*nz);
nx = nx/nr;
ny = ny/nr;
nz = nz/nr;
