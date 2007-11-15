function slices = slice_number2(lats,lons,latr,lonr,nproc,narc,lat0,nchunks,xi_width,eta_width,clon,clat,grot)
% this function get the slice number for a width 'lat0' belt between
% source and receiver points in the global mesh
% input: lats, lons, latr, lonr, nproc, narc
%        narc = 0 : minor arc
%        narc = 1 : major arc
%        lat0 : the width of the belt
% for 1-chunk or 2-chunk mesh, also need input:
%        nchunks, xi_width, eta_width, clon, clat, grot
%        
% output:  slice numbers 


% 
%clear all;
%lats = -13.82; lons=-67.25;
%latr = 18.79; lonr = 98.98;
%nproc = 5;

if (nargin ~= 7 && nargin ~= 13)
  disp('Number of arguments for slice_number should be either 6 or 12');
  return
end

if (nargin == 13 && nchunks > 1 && nproc ~= nproc_eta)
  disp('Number of processors in xi and eta directions need to be the same for nchunks > 1');
  return
end 

ths = (90-lats)/180.*pi; phs = lons/180.*pi;
thr = (90-latr)/180.*pi; phr = lonr/180.*pi;
xi_width = xi_width/180.*pi;
eta_width = eta_width/180.*pi;

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
      ph = linspace(ph1,ph2,Np);
  else
      ph = linspace(ph2-2*pi,ph1,Np);
  end
elseif (narc == 1) 
  if (delta < pi) 
      ph = linspace(ph2-2*pi,ph1,Np);
  else
      ph = linspace(ph1,ph2,Np);
  end
end

% calculate the rotation matrix from CHUNK_AB frame to the actual frame
if (nargin == 13); rotation_matrix = rotmat(clon,clat,grot); end
if (nargin == 7); nchunks = 6; xi_width = pi/2; yi_width = pi/2; end

% initialize
chunk(1:Np) = 10;
long_slice=[];

Nt = 11;
th0 = lat0 * pi/180;
th = linspace(-th0,th0,Nt);

for ii = 1 : Nt
for i = 1 : Np
    [thp(i),php(i)] = norm_rot_back(thn,phn,pi/2+th(ii),ph(i));
    [x,y,z] = tp2xyz(thp(i),php(i));

% rotate the points between source and receiver to CHUNK_AB standard coordinates
    if (nargin == 13)
      temp = [x,y,z] * rotation_matrix; 
      x=temp(1); y=temp(2); z=temp(3);
    end
    
    for k = 1: 6
        [xik,etak] = chunk_map(k,x,y,z);
        %[i,k,xik,etak]
         if (abs(xik) <= pi/4 && abs(etak) <= pi/4)
             chunk(i) = k;
             xi(i) = xik; eta(i) = etak;
         end
    end
% check chunk number
    if nargin == 12 && (chunk(i) > nchunks || abs(xi(i)) > xi_width/2 || abs(eta(i)) > eta_width/2)
      continue;
    end
    
    xi1 = xi(i) / xi_width * 2; eta1 = eta(i) / eta_width * 2;
    nproc_xi(i) = floor((xi1+1)/2 * nproc);
    nproc_eta(i) = floor((eta1+1)/2 * nproc);
    slice(i) = nproc * nproc * (chunk(i) - 1) + nproc * nproc_eta(i) + nproc_xi(i);
    info(i,:) = [chunk(i),nproc_xi(i),nproc_eta(i),slice(i),thp(i)*180/pi,php(i)*180/pi,xi(i),eta(i)];
    long_slice = [long_slice;slice(i)];
end
end

slices = compact_array(long_slice);

slices=slices';


    
 



