#!/usr/bin/perl

#
#  Script to create plate and continent boundaries in AVS format
#
#  Author : Dimitri Komatitsch, Caltech, August 2000
#

#
# first file contains continent boundaries with longitude and then latitude
#

# file with continents
  $name = 'continent_boundaries_gmt.dat';

  print STDOUT "Getting boundaries from file $name ...\n";

  open(FILEAVS,">AVS_continent_boundaries.inp");

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

  $pi = 3.14159265 ;

# write header for AVS (with point data)
      print FILEAVS "$numpoints $num_individual_lines 1 0 0\n";

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
      $theta = $pi/2. - atan2(0.99329534 * tan($latitude * $pi / 180.),1) ;
      $phi = $longitude * $pi / 180. ;

# compute the Cartesian position of the receiver (ignore ellipticity for AVS)
# assume a sphere of radius one
      $r_target = 1. ;
## DK DK make the radius a little bit bigger to make sure it is
## DK DK correctly superimposed to the mesh in final AVS figure
      $r_target = 1.015 ;
      $x_target = $r_target*sin($theta)*cos($phi) ;
      $y_target = $r_target*sin($theta)*sin($phi) ;
      $z_target = $r_target*cos($theta) ;

      print FILEAVS "$currentpoint $x_target $y_target $z_target\n";

      }
  }
  close(FILEGMT);

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
        print FILEAVS "$currentelem $currentline line $previouspoint $currentpoint\n";
        }
        $previous_was_comment = 0;
    }
  }
  close(FILEGMT);

  print FILEAVS " 1 1\n";
  print FILEAVS " Zcoord, meters\n";

# create data values for the points
  $currentpoint = 1;
  while($currentpoint <= $numpoints) {
      print FILEAVS "$currentpoint 255.\n";
      $currentpoint ++ ;
      }

  close(FILEAVS);


##
## DK DK define tangent function which is not standard in Perl
##
 sub tan { sin($_[0]) / cos($_[0])  }

