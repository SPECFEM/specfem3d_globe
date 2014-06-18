ETOPO2 Global 2-Minute Gridded Elevation Data Volume E1, 8/30/2001

Description:
      Topography for the entire world, including oceans, continents, and polar
      regions at a 2 minute resolution (~4 km) in meters.

GMT:
      This data set is included in the grdraster dataset for ease of plotting

File Format:
     The raw data is in ETOPO2v2c_i2_LSB.bin
     The file is a raw header-less binary stored as
     2-byte (16 bit) signed integers
     10800 columns and 5400 rows


Raw to Gridded Data:
    get file ETOPO2v2c_i2_LSB.zip from http://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2/ETOPO2v2-2006/ETOPO2v2c/raw_binary/
    unzip ETOPO2v2c_i2_LSB.zip
    xyz2grd ETOPO2v2c_i2_LSB.bin -Rd -I2m -Getopo2v2c.grd -F -ZTLh -V
    grd2xyz -Rg etopo2v2c.grd | awk '{ print $3 }' > topo_bathy_etopo2v2c_original_unmodified_unsmoothed.dat

Data Sources:
     For information about where the data originated check out
     et2src.htm and startet2.htm


-----------------------------------------------------------------------


This CD-ROM contains several elements: Software, Data, and Images:

GEODAS (GEOphysical DAta System) is an interactive database management
system developed by the National Geophysical Data Center (NGDC) for use
in the assimilation, storage and retrieval of geophysical data. The
GEODAS software is being used with several types of data including
marine trackline geophysical data, hydrographic (bathymetric) survey
data, aeromagnetic survey data, multibeam bathymetric data, and coastal
relief model bathymetric gridded data.

The GEODAS software which comes with this CD includes:
   Browser access to grids -- open startet2.htm in your browser after
       installing the GEODAS software from the "setup" directory
   NGDC Grid Translator Program
   Hydro-Plot Grid viewer Program

Data:

   Raw 2' grids in big-endian (ETOPO2.RAW) and little-endian (ETOPO2.dos)
        formats
        (16-bit signed integers, 10800 columns x 5400 rows)
   See Major Sources of Data" on startet2.htm for details
   For Geographic Information System (GIS) users, see the "ArcView Grids" link
        in the GEODAS Help Index

Images:

   4000x2000 color shaded relief JPEG image (grdet2.jpg)
   and 45-degree shaded relief image tiles with 2' or 5' pixels (1350x1350 pixels
   or 512x512 pixels), accessed through the STARTET2.HTM browser page (any browser
   and computer type)


********************
*** INSTALLATION ***
********************

MS Windows 95/98/NT:
To install see /setup/windows/readme.txt

UNIX Xwindows:
To install see /setup/xwindows/readme.txt

Macintosh users must have a PC emulation program to use the GEODAS features;
   the image display features of the STARTET2.HTM page work on any platform.

****************************
*** Technical Assistance ***
****************************

Dan Metzger
Dan.R.Metzger@noaa.gov
303-497-6542

Peter Sloss
Peter.W.Sloss@noaa.gov
303-497-6119

Dave Divins
David.Divins@noaa.gov
303-497-6505


Get the latest version of the software at our web site:
http://www.ngdc.noaa.gov/mgg/gdas/gx_announce.Html










