%% script by Hejun Zhu
%%
%% plots image of ETOPO1.xyz

clear all
close all

%% loads etopo1 array file
tmp = load('ETOPO1.xyz');

nx=21600;
ny=10800;

%% plots image
etopo=reshape(tmp,nx,ny);
pcolor(etopo(1:50:nx,1:50:ny))
shading interp
axis image
