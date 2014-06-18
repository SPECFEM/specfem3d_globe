%% script by Hejun Zhu
%%
%% converts ETOPO1_Ice_c.tif file to ETOPO1.xyz

clear all
close all


%% loads tiff file
fnm='ETOPO1_Ice_c.tif';
tmp = double(geotiffread(fnm));


[leny,lenx]=size(tmp);

%% prints out etopo1 array file
fid = fopen('ETOPO1.xyz','w');
for j = 1:leny
	for i = 1:lenx
		if i < lenx/2
			fprintf(fid,'%18.15f\n',tmp(j,lenx/2+i+1));
		else 
			fprintf(fid,'%18.15f\n',tmp(j,i-lenx/2+1));
		end 
	end 
end 
fclose(fid);

