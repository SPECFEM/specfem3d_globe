--------------------------------
readme
--------------------------------

VTK, The Visualization Toolkit
http://www.vtk.org/



- procedure to make a PNG image (flat earth):

------------------
NOTE: this requires VTK to be installed and compiled with python wrappers
      see small howto note below.
------------------


1. run specfem with the movie options (see Par_file):

   MOVIE_SURFACE = .true.
   MOVIE_COARSE  = .false.
   
   and adjust the time steps NSTEP_BETWEEN_FRAMES   
   
   this creates binary files in directory OUTPUT_FILES/ like: moviedata000100,...


2. convert binary files to GMT-files:

   in SPECFEM3D_GLOBE:  > make xcreate_movie_GMT_global
   
	         run it > ./xcreate_movie_GMT_global

                          choose option for binary (T)
                          to create individual files


3. create VTK files:

    > ./plot_movie_GMT_binary_VTK.pl OUTPUT_FILES/bin_movie_00***.d


4. render VTK files to create a PNG image:

    for each single file:
    
    > python plot_VTK.py OUTPUT_FILES/bin_movie_009000.d.vtk
    
    this creates a single PNG image 'bin_color.png'
    
    
optional, to add transparency:
    
    a) create a color image:

        > python plot_VTK.py OUTPUT_FILES/bin_movie_009000.d.vtk
        
       and a gray-scale image:        
       
        > python plot_VTK_gray.py OUTPUT_FILES/bin_movie_009000.d.vtk

    
    b) use the gray-scale image as alpha channel in the new file to
       create a file with transparency (opacity):

       > composite -compose CopyOpacity bin_mask.png bin_color.png bin_image.png
    
      this requires software installed from ImageMagick
      ( http://www.imagemagick.org/ )



------------------
VTK installation - python wrapper:
------------------

  download sources from:
  http://www.vtk.org/VTK/resources/software.html

  install with python wrappers and Geovis for example in /opt/vtk-5.4.2:
  > cd /opt/vtk-5.4.2
  > tar -xvf vtk-5.4.2.tar
  > cd VTK
  > ccmake .
  
    turn on options: BUILD_SHARED_LIBS ON
                     CMAKE_INSTALL_PREFIX /opt/vtk-5.4.2
                     VTK_USE_GEOVIS    ON
                     VTK_WRAP_PYTHON   ON
  > make
  > make install
   
  export your python path specifics for example in ~/.bashrc:
  
    # vtk python
    export PYTHONPATH=$PYTHONPATH:/opt/vtk-5.4.2/VTK/bin
    export PYTHONPATH=$PYTHONPATH:/opt/vtk-5.4.2/VTK/Wrapping/Python/
    export PATH=$PATH:/opt/vtk-5.4.2/VTK/bin    
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/vtk-5.4.2/VTK/lib/

  check if properly installed:
  > cd ~/SPECFEM3D_GLOBE/UTILS/VTK
  > python
  
    >>> from vtk import *
    >>> gs = vtkGeoProjection()
  
    if any of this fails, check your path default settings:
  
    >>> import sys
    >>> sys.path

    and fix the paths in your ~/.bashrc
  
  