function [thn,phn] = tp2norm(ths,phs,thr,phr)
% find out the normal of the plane constructed
% by (ths,phs), (thr,phr), and the origin using
% cross-product.

[xs,ys,zs] = tp2xyz(ths,phs);
[xr,yr,zr] = tp2xyz(thr,phr);
nx = ys*zr - zs*yr;
ny = -(xs*zr - zs*xr);
nz = xs*yr - ys*xr;

[thn,phn] = xyz2tp(nx,ny,nz);

