function [x,y,z] = xsection_translate(lats,lons,latr,lonr,narc,scale)
% this function figures out the x,y,z values for a translation within
% the source receiver cross-section with certain scale

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

%phr_new = phs_new + pi - 0.001;
pharray= sort([phs_new,phr_new]);
ph1 = pharray(1); ph2 = pharray(2);
delta = ph2-ph1;
%disp(strcat('Delta = ',num2str(delta * 180/pi)));

% we are looking at minor arc
if (narc == 0) 
  if (delta < pi)
      ph = (ph1+ph2+pi)/2;
  else
      ph = (ph2-2*pi+ph2-pi)/2;
  end
elseif (narc == 1) 
  if (delta < pi) 
      ph = (ph2-2*pi+ph2-pi)/2;
  else
      ph = (ph1+ph1+pi)/2;
  end
end

[thp,php] = norm_rot_back(thn,phn,pi/2.,ph);
[x,y,z] = tp2xyz(thp,php);
x = x*scale;
y = y*scale;
z = z*scale;
