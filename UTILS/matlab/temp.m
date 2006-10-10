% this function get the slice number for a width 'lat0' belt between
% source and receiver points in the 6-chunk global mesh
% input: lats, lons, latr, lonr, nproc, narc, lat0
%        narc = 0 : minor arc
%        narc = 1 : major arc
% output:  slice numbers

Nt = 11;
th0 = lat0 * pi/180;
th = linspace(-th0,th0,Nt);

for ii = 1 : Nt
for i = 1 : Np
    [thp(i),php(i)] = norm_rot_back(thn,phn,pi/2+th(ii).,ph(i));
