----------------------------------------------------------------------
PPM - point-profile models
----------------------------------------------------------------------

 for generic Vs models given as depth profiles at lon/lat using a text-file format like:

 #lon(deg), lat(deg), depth(km), Vs-perturbation wrt PREM(%), Vs-PREM (km/s)
  -10.00000       31.00000       40.00000      -1.775005       4.400000    
  -10.00000       32.00000       40.00000      -1.056823       4.400000    
 ...

   (the first line in the file is ignored)
   

1. in order to use it, create a symbolic link 'model.txt' 
   in this directory DATA/PPM/:
   > cd DATA/PPM/
   > ln -s my_ppm_model model.txt


2. change in DATA/Par_file the corresponding line to:
   #..
   MODEL                           = PPM
   .. 

3. any specifics can be added to the routine implementations in
   the subroutine file: 'model_ppm.f90'
   
   note: by default, the routine assumes that data given are Vs perturbations 
         and that the data have a regular spacing, i.e. constant increments, 
         for each lon,lat and depth column