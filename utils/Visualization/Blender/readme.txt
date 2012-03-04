--------------------------------
readme
--------------------------------

Blender, open source, cross platform suite of tools for 3D creation
http://www.blender.org/



- load example file 'bin_sphere.blend':

  run Blender:
  
      -> Menu -> File -> Open...
      
         select 'bin_sphere.blend'


  to render image: 
      
      -> Menu -> Render -> Render Current Frame


  it needs: ./topo_globe.jpg
            ./bin_movie_009000.d.png
               
  these images are used as textures on a sphere:  
  
  - a wavefield PNG image created with GMT:  
      
      > ../GMT/plot_movie_GMT_binary.pl OUTPUT_FILES/bin_movie_*.d
      
  - a topography image 'topo_globe.jpg':
      
      to create a topographic global map with GMT,

      1. create a grid file

        > grdraster ETOPO2 -Rg -I2m -Getopo2.grd

      2. create a gradient file

        > grdgradient etopo2.grd -Nt1 -A45 -Getopo2gradient.grd -V

      3. render image with gradient illumination
         (must copy topo_gray.cpt to current directory)

         projection uses: plate caree projection -JQ

        > grdimage etopo2.grd -Ietopo2gradient.grd -JQ0/0/15 -Ctopo_gray.cpt -V -P > topo_globe.ps

      4. convert to jpg using convert command from ImageMagick (http://www.imagemagick.org)
        
        > convert topo_globe.ps topo_globe.jpg

