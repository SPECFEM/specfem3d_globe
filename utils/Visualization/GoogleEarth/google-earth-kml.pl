#!/usr/bin/perl

use Getopt::Std;
use POSIX;

sub Usage {
  print STDERR <<EOF;

Usage: google-earth-kml.pl file_name


  ex. ./google-earth-kml.pl gmt_movie*.gif

EOF
exit(1)
}


if (!getopts('l:L:cdts')) {die('Check input arguments\n');}
@ARGV > 0 or Usage();

#######################################################
## USER/EVENT PARAMETERS

# look at viewpoint
$lookat_lon= 78.0;
$lookat_lat= 10.0;

# from CMTSOLUTION file
$starttime="2003-12-26T01:56:52.4Z";

# time increment per image ( 100 time steps per 0.124 s ~ 12.4 s)
$tinc = 12.4;

#######################################################

open(KML,">animation.kml") or die("error file open");

print KML '<?xml version="1.0" encoding="UTF-8"?>';print KML "\n";
print KML '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">';print KML "\n";
print KML "<Folder>\n";
print KML "<name>Animation</name>\n";
print KML "<open>0</open>\n";


# look at view
print KML "<LookAt>\n";
print KML "<longitude> $lookat_lon </longitude>\n";
print KML "<latitude> $lookat_lat </latitude>\n";
print KML "<altitude> 0 </altitude>\n";
print KML "<range>10929687.63258196</range>\n";
print KML "<tilt>0</tilt>\n";
print KML "<heading>6.327702684212388</heading>\n";
print KML "<altitudeMode>relativeToGround</altitudeMode>\n";
print KML "</LookAt>\n";


foreach $file (@ARGV) {
  print "processing $file\n";
  if (! -f $file) {die (" check to see if $file exists or not\n");}

  # image overlay
  print KML "<GroundOverlay>\n";
  print KML " <name>test-movie</name>\n";
  print KML " <TimeSpan>\n";


  # takes starttime from last image
  ($cal,$utc) = split("T",$starttime);
  ($year,$month,$day) = split("-",$cal);

  $utc2 = substr($utc,1,10);
  ($hh,$mm,$ss) = split(":",$utc2);

  # adds time increment
  $ss = $ss+$tinc;

  # calculates image time in hours/minutes/seconds
  if( $ss > 60.0 ){
    $ss = $ss - 60.0;
    $mm = $mm + 1;
  }
  if( $mm > 60 ){
    $hh = $hh + 1;
    $mm = $mm - 60;
  }
  if( $hh > 23 ){
    $hh = 00;
    $day = $day+1;
  }

  $time = sprintf("%4.4i-%2.2i-%2.2iT%2.2i:%2.2i:%03.1fZ",$year,$month,$day,$hh,$mm,$ss);

  print KML " <begin>$starttime</begin>\n";
  print KML " <end>$time</end>\n";
  print KML " </TimeSpan>\n";
  print KML " <color>a1ffffff</color>\n";
  print KML " <Icon>\n";
  print KML "  <href>$file</href>\n";
  print KML "  <viewBoundScale>0.75</viewBoundScale>\n";
  print KML " </Icon>\n";

  # image box (spans the whole globe)
  print KML " <LatLonBox>\n";
  print KML "  <north>90.0</north>\n";
  print KML "  <south>-90.0</south>\n";
  print KML "  <east>180.0</east>\n";
  print KML "  <west>-180.0</west>\n";

  print KML " </LatLonBox>\n";
  print KML "  </GroundOverlay>\n";

  # sets new starttime
  $starttime = $time;

}

print KML "</Folder>\n";
print KML "</kml>";


close(KML);

print "   Done !\n";
