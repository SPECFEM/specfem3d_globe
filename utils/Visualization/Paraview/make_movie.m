%
% make_movie.m
% /SPECFEM3D_GLOBE/UTILS/Visualization/Paraview/make_movie.m
%
% template script for using matlab to make a movie from a set of jpgs 
%
% NOTE: VideoWriter requires version R2010a or later (type "version")
%
% Alternatively, try the command-line using ImageMagick:
%    convert -delay 8 -resize 400x320 *.jpg movie.mpg
%
%

close
clear all

%----------------
% USER PARAMETERS

ifiles = '*.jpg';           % input files
ofile = 'simulation.avi';   % output file
fps = 2;                    % frames per second

% optional: construct text label for each jpg
itextlabel = 1;             % =1 to add text label (see below)
frames = [1000:400:17800]'; % frames
DT = 7.6219589E-02;         % time step
t0 = -60;                   % origin time (USER_T0)
tvec = (t0 + DT*frames)/60;
tx0 = 30; ty0 = 40;         % position of a text label (in pixels)

%----------------

list = dir(ifiles);         % create structure
files = {list.name};

vid_out = VideoWriter(ofile);
vid_out.FrameRate = fps;
open(vid_out);

fig = figure(1); clf

% loop over each input file
for i = 1:length(files)
    img = imread(files{i});
    imshow(img, 'Border', 'tight', 'InitialMagnification', 80);
    
    % add text label
    if itextlabel==1
        stlab = sprintf('Time = %3.0f minutes',tvec(i));
        text(tx0,ty0,stlab,'fontsize',18,'color','w','fontweight','bold');
    end
    
    drawnow
    writeVideo(vid_out, getframe(fig));
end

close(vid_out);

%=========================================================
