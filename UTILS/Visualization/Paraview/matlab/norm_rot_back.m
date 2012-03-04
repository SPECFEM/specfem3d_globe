function [th_old,ph_old] = norm_rot_back(thn,phn,th,ph)
% coordinates change from (th_old,ph_old) to (th,ph) 
% according to a rotation that converts (thn,phn) to
% z axis

rot=[cos(thn)*cos(phn),cos(thn)*sin(phn),-sin(thn);
    -sin(phn),cos(phn),0;
    sin(thn)*cos(phn),sin(thn)*sin(phn),cos(thn)];
[x,y,z] = tp2xyz(th,ph);
temp=rot' * [x,y,z]';
x_old = temp(1); y_old = temp(2); z_old = temp(3);
[th_old,ph_old] = xyz2tp(x_old,y_old,z_old);

