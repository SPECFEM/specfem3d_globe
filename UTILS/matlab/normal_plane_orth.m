function [xc,yc,zc,nx,ny,nz] = normal_plane_orth(lats,lons,latr,lonr)
% find out the coord (xc,yc,zc) of the center of the vector source-receiver
% and the normalised coordinates (nx,ny,nz) of this vector 
% for a "doughnut cross-section"

ts = (90 - lats) * pi/180;
ps = lons * pi / 180;
tr = (90 - latr) * pi/180;
pr = lonr * pi / 180;

[xs,ys,zs] = tp2xyz(ts,ps);
[xr,yr,zr] = tp2xyz(tr,pr);

nx=xr-xs;
ny=yr-ys;
nz=zr-zs;
nr=sqrt(nx*nx+ny*ny+nz*nz);
nx = nx/nr;
ny = ny/nr;
nz = nz/nr;

xc=(xr+xs)/2;
yc=(yr+ys)/2;
zc=(zr+zs)/2;

