--------------------------------
readme
--------------------------------

GMT, the Generic Mapping Tools 
http://www.soest.hawaii.edu/GMT/ 



- procedure to make a GMT image:

1. run specfem with the movie options (see Par_file):

   MOVIE_SURFACE = .true.
   MOVIE_COARSE  = .false.
   
   and adjust the time steps NSTEP_BETWEEN_FRAMES   
   
   this creates binary files in directory OUTPUT_FILES/ like: moviedata000100,...


now, you have basically two options:


2.A) convert binary files to GMT-files:

   in SPECFEM3D_GLOBE:  > make xcreate_movie_AVS_DX
                          
	         run it > ./xcreate_movie_AVS_DX

                          choose option 4 for individual files
                          
   outputs file like: OUTPUT_FILES/gmt_movie_000100.xyz
   

  - render gif files using GMT:

    > ./plot_gmt_movie_file.pl OUTPUT_FILES/gmt_movie*.xyz



2.B) convert binary files to GMT-files:

   in SPECFEM3D_GLOBE:  > make xcreate_movie_GMT_global
   
	         run it > ./xcreate_movie_GMT_global

                          choose either option for ascii (F) or binary (T)
                          to create individual files


  - render ascii files to create ps files using GMT:

    > ./plot_movie_GMT_ascii.pl OUTPUT_FILES/ascii_movie_00***.d


  - render binary files to create PNG files using GMT:

    > ./plot_movie_GMT_binary.pl OUTPUT_FILES/bin_movie_00***.d


please, change the rendering scripts plot_movie_GMT_ascii.pl, etc.
according to your needs.

note: these scripts make use of commands from ImageMagick (http://www.imagemagick.org)