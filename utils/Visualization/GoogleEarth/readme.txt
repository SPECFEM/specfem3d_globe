--------------------------------
readme
--------------------------------

Google Earth
http://earth.google.com/



- procedure to make a google earth movie plot:

1. run specfem with the movie options (see Par_file):

   MOVIE_SURFACE = .true.
   MOVIE_COARSE  = .false.
   
   and adjust the time steps NSTEP_BETWEEN_FRAMES   
   
   this creates binary files in directory OUTPUT_FILES/ like: moviedata000100,...

   
2. convert binary files to GMT-files:

   in SPECFEM3D_GLOBE:  > make xcreate_movie_AVS_DX
                          
	         run it > ./xcreate_movie_AVS_DX

                          choose option 4 for individual files
                          
   outputs file like: OUTPUT_FILES/gmt_movie_000100.xyz
   

3. render gif files using GMT:

    > UTILS/Visualization/GMT/plot_gmt_movie_file.pl OUTPUT_FILES/gmt_movie*.xyz


4. create google animation.kml file:

   > ./google-earth-kml.pl OUTPUT_FILES/gmt_movie*.gif
     
     
     this will create a file 'animation.kml' which can be opened
     with GoogleEarth application (http://earth.google.com/)


please, change the rendering script 'google-earth-kml.pl' according to your needs.     
