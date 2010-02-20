function [nx,ny,nz] = tp2norm(ths,phs,thr,phr)
% find out the normal of the plane (nx,ny,nz) 
% constructed by (ths,phs), (thr,phr), and the origin.

[xs,ys,zs] = tp2xyz(ths,phs);
[xr,yr,zr] = tp2xyz(thr,phr);
nx = ys*zr - zs*yr;
ny = -(xs*zr - zs*xr);
nz = xs*yr - ys*xr;
nr = sqrt(nx*nx+ny*ny+nz*nz);
nx = nx/nr;
ny = ny/nr;
nz = nz/nr;

