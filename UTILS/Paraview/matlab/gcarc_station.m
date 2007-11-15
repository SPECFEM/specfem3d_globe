function [thp,php] = gcarc_station(lats,lons,latr,lonr,narc)
% this function get the geographic coordinates for points along
% the great circle path between the
% source and receiver points in the 6-chunk global mesh
% input: lats, lons, latr, lonr, nproc, narc
%        narc = 0 : minor arc
%        narc = 1 : major arc
% output:  thp(1:180), php(1:180), 1 degree apart


% 
%clear all;
%lats = -13.82; lons=-67.25;
%latr = 18.79; lonr = 98.98;
%nproc = 5;


ths = (90-lats)/180.*pi; phs = lons/180.*pi;
thr = (90-latr)/180.*pi; phr = lonr/180.*pi;

[thn,phn] = tp2norm(ths,phs,thr,phr);
[ths_new,phs_new] = norm_rot(thn,phn,ths,phs);
[thr_new,phr_new] = norm_rot(thn,phn,thr,phr);

if (ths_new-pi/2) > 1e-2
    disp('new lat of source is not 0'); return;
end
if (thr_new-pi/2) > 1e-2
    disp('new lat of receiver is not 0');return;
end

Np = 180;
%phr_new = phs_new + pi - 0.001;
pharray= sort([phs_new,phr_new]);
ph1 = pharray(1); ph2 = pharray(2);
delta = ph2-ph1;
%disp(strcat('Delta = ',num2str(delta * 180/pi)));

% we are looking at minor arc
if (narc == 0) 
  if (delta < pi)
      ph = linspace(ph1,ph2+pi,Np);
  else
      ph = linspace(ph2-2*pi,ph2-pi,Np);
  end
elseif (narc == 1) 
  if (delta < pi) 
      ph = linspace(ph2-2*pi,ph2-pi,Np);
  else
      ph = linspace(ph1,ph1+pi,Np);
  end
end

for i = 1 : Np
    [thp(i),php(i)] = norm_rot_back(thn,phn,pi/2.,ph(i));
end
