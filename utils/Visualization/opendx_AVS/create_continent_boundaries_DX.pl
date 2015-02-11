#!/usr/bin/perl

#
#  Script to create plate and continent boundaries in OpenDX format
#
#  Author : Dimitri Komatitsch, Caltech, August 2000
#

#
# first file contains continent boundaries with longitude and then latitude
#

# file with continents
  $name = 'continent_boundaries_gmt.dat';

  print STDOUT "Getting boundaries from file $name ...\n";

  open(FILEAVS,">DX_continent_boundaries.dx");

  $numsegments = `cat $name | grep ">" | wc -l `;
  chop $numsegments;
  print STDOUT "There are $numsegments contours\n";

  $numpoints = `cat $name | grep -v ">" | wc -l `;
  chop $numpoints;
  print STDOUT "There are $numpoints data points\n";

# read the GMT file to get the number of individual line segments
  $currentelem = 0;
  $previous_was_comment = 1;
  open(FILEGMT,"<$name");
  while($line = <FILEGMT>) {
    chop $line;
#   get line marker (comment in file)
    if(substr($line,0,1) eq '>') {
      $previous_was_comment = 1;
      }
    else {
      if($previous_was_comment == 0) {
        $currentelem ++;
        }
        $previous_was_comment = 0;
    }
  }
  close(FILEGMT);
  $num_individual_lines = $currentelem;
  print STDOUT "There are $num_individual_lines individual line segments\n";

  $pi = 3.1415926535 ;

# write header for AVS (with point data)
      print FILEAVS "object 1 class array type float rank 1 shape 3 items $numpoints data follows\n";

# read the GMT file to get the points
  $currentpoint = 0;
  open(FILEGMT,"<$name");
  while($line = <FILEGMT>) {
    chop $line;
#   get point only if line is not a comment
    if(substr($line,0,1) ne '>') {
      $currentpoint ++;

# longitude is the number before the white space
 $longitude = substr($line,0,index($line," "));

# latitude is the number after the white space
 $latitude = substr($line,index($line," ")+1);

# convert geographic latitude to geocentric colatitude and convert to radians
#     $theta = $pi/2. - atan2(0.99329534 * tan($latitude * $pi / 180.),1) ;
      $theta = $pi/2. - $latitude * $pi / 180. ;
      $phi = $longitude * $pi / 180. ;

# compute the Cartesian position of the receiver (ignore ellipticity)
# assume a sphere of radius one
## DK DK make the radius a little bit bigger to make sure it is
## DK DK correctly superimposed to the mesh in final AVS or OpenDX figure
#     $r_target = 1.015 ;
      $r_target = 1.007 ;
      $x_target = $r_target*sin($theta)*cos($phi) ;
      $y_target = $r_target*sin($theta)*sin($phi) ;
      $z_target = $r_target*cos($theta) ;

      print FILEAVS "$x_target $y_target $z_target\n";

      }
  }
  close(FILEGMT);

 print FILEAVS "object 2 class array type int rank 1 shape 2 items $num_individual_lines data follows\n";

# read the GMT file to get the lines
  $currentline = 0;
  $currentelem = 0;
  $currentpoint = 0;
  $previous_was_comment = 1;
  open(FILEGMT,"<$name");
  while($line = <FILEGMT>) {
    chop $line;
#   get line marker (comment in file)
    if(substr($line,0,1) eq '>') {
      $currentline ++;
      $currentpoint ++;
      $previous_was_comment = 1;
      print STDOUT "processing contour $currentline named $line\n";
      }
    else {
      if($previous_was_comment == 0) {
        $previouspoint = $currentpoint;
        $currentelem ++;
        $currentpoint ++;
        $previouspointdx = $previouspoint - 1;
        $currentpointdx = $currentpoint - 1;
        print FILEAVS "$previouspointdx $currentpointdx\n";
        }
        $previous_was_comment = 0;
    }
  }
  close(FILEGMT);

 print FILEAVS "attribute \"element type\" string \"lines\"\n";
 print FILEAVS "attribute \"ref\" string \"positions\"\n";
 print FILEAVS "object 3 class array type float rank 0 items $numpoints data follows\n";

# create data values for the points
  $currentpoint = 1;
  while($currentpoint <= $numpoints) {
      print FILEAVS "255.\n";
      $currentpoint ++ ;
      }

 print FILEAVS "attribute \"dep\" string \"positions\"\n";
 print FILEAVS "object \"irregular connections  irregular positions\" class field\n";
 print FILEAVS "component \"positions\" value 1\n";
 print FILEAVS "component \"connections\" value 2\n";
 print FILEAVS "component \"data\" value 3\n";
 print FILEAVS "end\n";

  close(FILEAVS);

##
## DK DK define tangent function which is not standard in Perl
##
 sub tan { sin($_[0]) / cos($_[0])  }

