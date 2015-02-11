#!/usr/bin/perl
#
# reads in ascii files produced from xcreate_movie_GMT_global
# and creates PNG files with transparency
#

use Getopt::Std;
use POSIX;

sub Usage {
  print STDERR <<EOF;

Usage: plot_movie_GMT_ascii_PNG.pl file_name

  ex. ./plot_movie_GMT_ascii_PNG.pl OUTPUT_FILES/ascii_movie_*.d

  converts ps-file to *.png file format with transparency

EOF
exit(1)
}

if (!getopts('l:L:cdts')) {die('Check input arguments\n');}
@ARGV > 0 or Usage();


#######################################################
## PARAMETERS

# point locations
$xy_file = "OUTPUT_FILES/ascii_movie.xy";

# global region
$R = "-Rd";

#plate carre projection
$JM = "-JQ0/0/15";

# interpolation on a regular spaced grid file
$mesh_degree = "0.8/0.8";
$intp_degree = "0.1/0.1";

# verbosity
$verbose = ""; # "-V";
# grid
$grdfile = "ascii_movie.grd";
# annotation
$B = " -B10/10wesn ";

# color tables
$colorcpt = "blue_white_red.cpt";
$maskcpt = "gray_pyramid_inv.cpt";

#######################################################

open FILE, "$xy_file" or die $!;
@lines = <FILE>;
close(FILE);
$nlines = @lines;


open(CSH,">plot_movie.csh");
print CSH "gmtset BASEMAP_TYPE plain ANOT_FONT_SIZE 9 HEADER_FONT_SIZE 15\n";


foreach $file (@ARGV) {

  if (not -f $file) {die("No $file\n");}

  print "Processing frame $file...\n";

  # reads displacement file
  open FILE, "$file" or die $!;
  @lines_f = <FILE>;
  close(FILE);
  $nlines_f = @lines_f;

  if ($nlines_f != $nlines) {die("number of lines differ\n");}

  # determines min/max of displacements
  $minmax = `minmax $file | awk '{print \$5}' | sed -e "s:<: :" | sed -e "s:>: :" `;
  chomp($minmax);
  ($min,$max) = split("/",$minmax);

  # transforms displacements to range [0,255] for colors
  if( abs($min) > abs($max) ){$max = abs($min);}
  open FILE, ">$file.xyz" or die $!;
  for($i=0;$i<$nlines;$i++){
    # scale between 0,1
    $val = ($lines_f[$i] + $max)/(2.0*$max);
    # scale between 0, 255
    $val = $val * 255.0;

    # prints value together with coordinates
    $coord = $lines[$i];
    chomp($coord);
    #if( $val < 107 || $val > 148) {print FILE "$coord $val \n";}
    print FILE "$coord $val \n";
  }
  close(FILE);

  # interpolates values on a regular spaced grid file
  print CSH "xyz2grd $file.xyz -I$mesh_degree $R -G$grdfile  $verbose \n"; # -N127.5
  print CSH "grdsample $grdfile  -I$intp_degree -G$grdfile.1 -F $verbose\n";


  $ps_file = "$file.ps";

  # start ps-file
  print CSH "psxy -JX1/1 -R0/1/0/1 -K -P $verbose -X1 -Y1 <<EOF >$ps_file\nEOF\n";
  print CSH "grdimage $grdfile.1 $JM $R -Q -C$colorcpt $B -K -O $verbose -P >> $ps_file\n";
  print CSH "psxy -JX1/1 -R0/1/0/1 -O -P $verbose <<EOF >>$ps_file\nEOF\n";
  # end ps-file

  # start ps-file mask
  print CSH "psxy -JX1/1 -R0/1/0/1 -K -P $verbose -X1 -Y1 <<EOF >$ps_file.mask\nEOF\n";
  print CSH "grdimage $grdfile.1 $JM $R -Q -C$maskcpt $B -K -O $verbose -P >> $ps_file.mask\n";
  print CSH "psxy -JX1/1 -R0/1/0/1 -O -P $verbose <<EOF >>$ps_file.mask\nEOF\n";
  # end ps-file mask

  # convert to png
  print CSH "convert -crop 1080x540+6+6 $ps_file $file.disp.png\n";
  print CSH "convert -crop 1080x540+6+6 $ps_file.mask $file.mask.png\n";

  # creates file with transparency (opacity)
  print CSH "composite -compose CopyOpacity $file.mask.png $file.disp.png $file.png\n";

  print CSH "\necho $file.png\n";

  # remove temporary data file
  print CSH "rm -f $file.xyz\n";

  # removes temporary images
  print CSH "rm -f $ps_file $ps_file.mask $file.mask.png $file.disp.png \n";

}
close(CSH);

print "\nplotting... \n\n";

system("csh -fv plot_movie.csh");
system("rm -f plot_movie.csh $grdfile $grdfile.1");

