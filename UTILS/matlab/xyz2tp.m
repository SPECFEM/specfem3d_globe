function [th,ph] = xyz2tp(x,y,z)
% convert x,y,z to a point on unit sphere

ph = atan2(y,x);
th = atan2(sqrt(x*x+y*y),z);
if (th < 0) 
    th = th + pi;
end
if (ph < 0)
    ph = ph + 2*pi;
end
