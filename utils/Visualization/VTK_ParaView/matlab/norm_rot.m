function [th_new,ph_new] = norm_rot(thn,phn,th,ph)
% coordinates change from (th,ph) to (th_new,ph_new) 
% according to a rotation that converts (thn,phn) to
% z axis

rot=[cos(thn)*cos(phn),cos(thn)*sin(phn),-sin(thn);
    -sin(phn),cos(phn),0;
    sin(thn)*cos(phn),sin(thn)*sin(phn),cos(thn)];
[x,y,z] = tp2xyz(th,ph);
temp=rot * [x,y,z]';
x_new = temp(1); y_new = temp(2); z_new = temp(3);
[th_new,ph_new] = xyz2tp(x_new,y_new,z_new);


