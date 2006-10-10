function [x,y,z] = tp2xyz(th,ph)
% convert (th,ph) to (x,y,z) on unit sphere
x = sin(th) * cos(ph);
y = sin(th) * sin(ph);
z = cos(th);