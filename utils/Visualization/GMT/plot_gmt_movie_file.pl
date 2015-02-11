#!/usr/bin/perl
#
# reads in files produced from xcreate_movie_AVS_DX for GMT plotting
# and creates gif files (with option for one-color transparency)
#

use Getopt::Std;
use POSIX;

sub Usage {
  print STDERR <<EOF;

Usage: plot_gmt_movie_file.pl file_name

  ex. ./plot_gmt_movie_file.pl OUTPUT_FILES/gmt_movie_000***.xyz

EOF
exit(1)
}

if (!getopts('l:L:cdts')) {die('Check input arguments\n');}
@ARGV > 0 or Usage();


#######################################################
## region
# global
$R = "-Rd";

#plate carre projection (centered at 0 meridian, parallel to equator, plot width 15 inches)
$JM = "-JQ0/0/15";

# use white as transparent color
$transparency = 0;

#######################################################


open(CSH,">plot_movie.csh");
print CSH "gmtset BASEMAP_TYPE plain ANOT_FONT_SIZE 9 HEADER_FONT_SIZE 15\n";
print CSH "makecpt -Cpolar -T0/255/1 -Z -V > grd.cpt\n";

foreach $file (@ARGV) {

  if (not -f $file) {die("No $file\n");}

  print "Processing frame $file...\n";

  open FILE, "$file" or die $!;
  @lines_f = <FILE>;
  close(FILE);
  $nlines_f = @lines_f;

  $ps_file = "$file.ps";

  $grdfile = "ascii_movie.grd";
  $B = " -B10/10wesn ";

  # start ps-file
  print CSH "psxy -JX1/1 -R0/1/0/1 -K -P -V -X1 -Y1 <<EOF >$ps_file\nEOF\n";


  print CSH "xyz2grd $file -I0.4/0.4 $R -G$grdfile -N127.5 -V \n";
  print CSH "grdsample $grdfile -G$grdfile.1 -I0.4/0.4 -F -V\n";
  print CSH "grdimage $grdfile.1 $JM $R  -Cgrd.cpt $B -K -O -V -P >> $ps_file\n";

  # coast
####  print CSH "pscoast $JM $R -W0.1 -Dl  -A10000 -K -O -P -V >> $ps_file \n";
  print CSH "pscoast $JM $R -W0.1 -Dh -A1000 -K -O -P -V >> $ps_file \n";
  # color scale
  print CSH "psscale -D3/-0.5/3/0.2h -Ba50:'': -Cgrd.cpt -K -O -V -P >> $ps_file\n";
  # base map
  print CSH "psbasemap  -R -J -Ba90/a30:.'':WeSn -O -K -P -V  >> $ps_file\n";

  # end ps-file
  print CSH "psxy -JX1/1 -R0/1/0/1 -O -P -V <<EOF >>$ps_file\nEOF\n";

  # convert to gif
  if( $transparency == 1 ){
    print CSH "convert -crop 1080x540+72+6 -transparent 'rgb(254,254,255)' $ps_file $file.gif\n";
  }
  else{
    print CSH "convert -crop 1080x540+72+6 $ps_file $file.gif\n";
  }

  print CSH "rm -f $ps_file\n";

  print CSH "\necho $file.gif\n";
}
close(CSH);

print "\nplotting... \n\n";

system("csh -fv plot_movie.csh");
system("rm -f plot_movie.csh $grdfile $grdfile.1 grd.cpt");

