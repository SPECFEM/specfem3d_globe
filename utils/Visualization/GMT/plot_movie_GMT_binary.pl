#!/usr/bin/perl
#
# reads in binary files produced from xcreate_movie_GMT_global
# and creates ps files
#

use Getopt::Std;
use POSIX;

sub Usage {
  print STDERR <<EOF;

Usage: plot_movie_GMT_binary.pl file_name

  ex. ./plot_movie_GMT_binary.pl OUTPUT_FILES/bin_movie_*.d

EOF
exit(1)
}

if (!getopts('l:L:cdts')) {die('Check input arguments\n');}
@ARGV > 0 or Usage();


#######################################################
## PARAMETERS

# point locations
$xy_file = "OUTPUT_FILES/bin_movie.xy";


# global region
$R = "-Rd";

#plate carre projection
$JM = "-JQ0/0/15";

# interpolation
$interp = "-I1.5/1.5";
#######################################################

# verbosity
$verbose = ""; # "-V";
# grid
$grdfile = "movie.grd";
# annotation
$B = " -B10/10wesn ";

open(CSH,">plot_movie.csh");
print CSH "gmtset BASEMAP_TYPE plain ANOT_FONT_SIZE 9 HEADER_FONT_SIZE 15\n";
print CSH "makecpt -Cpolar -T0/255/2 -Z -V > grd.cpt\n";

print "Processing locations: $xy_file\n";

# reads in locations (binary data)
open FILE, "$xy_file" or die $!;
binmode FILE;

@lines= () ;
$line="";
$is_ok = 1;
$count = 0;
while ( $is_ok != 0) {
  read FILE, $junk, 4;

  $n1 = read FILE, $data1, 4;
  if( $n1 != 0 ){
    ($lon) = unpack( "f",$data1);
  }
  $n2 = read FILE, $data2, 4;
  if( $n2 != 0 ){
    ($lat) = unpack( "f",$data2);
  }
  read FILE, $junk, 4;

  #print "$n1 $n2 bytes read: lon=$lon \t \t lat=$lat\n";

  if( $n1 != 0 && $n2 != 0  ){
    $line = "$lon $lat \n";
    push(@lines, $line);
    $is_ok = 1;
  }
  else{ $is_ok = 0;}

  $count++;
  #if( $count == 8 ) {die("end");}
}
$nlines = @lines;
close(FILE);
print "  number of lines: $nlines\n\n";


foreach $file (@ARGV) {

  if (not -f $file) {die("No $file\n");}

  print "Processing frame $file...\n";

  # reads displacement file (binary data)
  open FILE, "$file" or die $!;
  binmode FILE;

  @lines_f= () ;
  $line="";
  $is_ok = 1;
  $min = 1.e30;
  $max = -1.e30;
  while ( $is_ok != 0) {
    read FILE, $junk, 4;

    $n1 = read FILE, $data1, 4;
    if( $n1 != 0 ){
      ($val) = unpack( "f",$data1);
    }

    read FILE, $junk, 4;

    if( $power_scaling > 0 ) {
      if( $val > 0 ){ $val = $val**$power_scaling;
      }else{
        $val = - abs($val)**$power_scaling;
      }

    }

    #print "$n1 bytes read: value=$val \n";

    if( $n1 != 0 ){
      $line = "$val \n";
      push(@lines_f, $line);

      # determines min/max of displacements
      if( $val < $min) {$min = $val;}
      if( $val > $max) {$max = $val;}

      $is_ok = 1;
    }
    else{ $is_ok = 0;}
  }
  $nlines_f = @lines_f;
  close(FILE);

  if ($nlines_f != $nlines) {die("number of lines differ\n");}


  # transforms displacements to range [0,255] for colors
  if( abs($min) > abs($max) ){$max = abs($min);}
  open FILE, ">$file.xyz" or die $!;

  @elem_corners = ();
  @elem_val = ();
  for($i=0;$i<$nlines;$i++){
    # scale between 0,1
    $val = ($lines_f[$i] + $max)/(2.0*$max);
    # scale between 0, 255
    $val = $val * 255.0;

    # prints value together with coordinates
    $coord = $lines[$i];
    chomp($coord);

    # adds point
    print FILE "$coord $val \n";

    # stores values in arrays
    push(@elem_val, $val );
    $line = "$coord $val ";
    push(@elem_corners, $line);
  }
  close(FILE);

  # output file
  $ps_file = "$file.ps";

  # interpolates displacement field
  print CSH "xyz2grd $file.xyz $interp $R -G$grdfile -N127.5 -V \n";
  print CSH "grdsample $grdfile -G$grdfile.1 $interp -F -V\n";
  print CSH "grdimage $grdfile.1 $JM $R  -Cgrd.cpt $B -K -V -P > $ps_file\n";

  # draws grid points
  print CSH "awk '{print \$0}' $file.xyz | psxy -J -R -Sc0.01 -Cgrd.cpt -V -K -O -P >> $ps_file \n";

  # river, states, coast
  #print CSH "pscoast $JM $R -W4 -Na -Dh -K -O -P -V >> $ps_file \n";

  # coast only
  print CSH "pscoast $JM $R -W1 -Dl  -A10000 -K -O -P -V >> $ps_file \n";

  # color scale
  print CSH "psscale -D3/-0.5/3/0.2h -Ba50:'': -Cgrd.cpt -K -O -V -P >> $ps_file\n";

  # base map
  print CSH "psbasemap  -R -J -Ba90/a30:.'':WeSn -O -P -V  >> $ps_file\n";

  # remove temporary data file
  print CSH "rm -f $file.xyz \n";

  # removes temporary images

  print CSH "echo \n";
  print CSH "echo 'plotted: $file.ps' \n";
  print CSH "echo \n";

}
close(CSH);

print "\nplotting... \n\n";

# executes script
system("csh -f plot_movie.csh");

# cleanup
system("rm -f plot_movie.csh");
system("rm -f $grdfile $grdfile.1 grd.cpt");
