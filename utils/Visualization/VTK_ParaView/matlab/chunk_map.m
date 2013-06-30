function [xi,eta] = chunk_map(k,xx,yy,zz)
% this program get the xi,eta for (xx,yy,zz)
% point under the k'th chunk coordinate
% transformation
small = 1d-6;
x = xx; y = yy; z = zz;
if (x < small & x >= 0); x = small; end
if (x > -small & x < 0); x = -small; end
if (y < small & y >= 0); y = small; end
if (y > -small & y < 0); y = -small; end
if (z < small & z >= 0); z = small; end
if (z > -small & z < 0); z = -small; end

if k == 1 % CHUNK_AB 
  xi = atan(y/z); eta = atan(-x/z);
  if (z < 0); xi = 10; end
elseif k == 2 % CHUNK_AC
  xi = atan(-z/y); eta = atan(x/y);
  if (y > 0); xi = 10; end
elseif k == 3 % CHUNK_BC
  xi = atan(-z/x); eta = atan(-y/x);
  if (x > 0); xi = 10; end
elseif k == 4 % CHUNK_AC'
  xi = atan(-z/y); eta = atan(-x/y);
  if (y < 0); xi = 10; end
elseif k == 5 % CHUNK_BC'
  xi = atan(z/x); eta = atan(-y/x);
  if (x < 0); xi = 10; end
elseif k == 6 % CHUNK_AB'
  xi = atan(y/z); eta = atan(x/z);
  if (z > 0); xi = 10; end
else
  disp('k != 1 - 6');
  return;
end

