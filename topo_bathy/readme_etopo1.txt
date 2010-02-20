----------------------------------------------------------------------
ETOPO1  
----------------------------------------------------------------------


1. download the whole world gridfile, ETOPO1 Ice Surface, in (georeferenced) tiff format 
   from this website:
    
   http://www.ngdc.noaa.gov/mgg/global/global.html

   and rename file to 'ETOPO1_Ice_c.tif'
   

2. process with matlab script 'xprocess_ETOPO1.m' to create file 'ETOPO1.xyz':

  > matlab -nojvm -r "xprocess_ETOPO1"

  note: this will probably take a while (~up to a few hours) and produce a rather big file (~4.7GB).


3. uncomment the lines in file 'constants.h' which refer to ETOPO1 as topography 


in order to check the newly created file 'ETOPO1.xyz', a plotting 
script 'xread_ETOPO1.m' is provided to produce an image plot of the file. 
  
  it can be called by:
  > matlab -nojvm -r "xread_ETOPO1"



ETOPO1 reference:

Amante, C. and B. W. Eakins,
ETOPO1 1 Arc-Minute Global Relief Model: 
Procedures, Data Sources and Analysis. 
NOAA Technical Memorandum NESDIS NGDC-24, 19 pp, March 2009.