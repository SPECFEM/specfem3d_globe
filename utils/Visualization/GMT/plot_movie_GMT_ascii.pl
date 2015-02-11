#!/usr/bin/perl
#
# reads in ascii files produced from xcreate_movie_GMT_global
# and creates ps files
#

use Getopt::Std;
use POSIX;

sub Usage {
  print STDERR <<EOF;

Usage: plot_movie_GMT_ascii.pl file_name

  ex. ./plot_movie_GMT_ascii.pl OUTPUT_FILES/ascii_movie_000100.d

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
$JM = "-JQ0/0/5";

#interpolation
$interp = "-I1.5/1.5";
#######################################################

# parameters
$grdfile = "movie.grd";
$B = " -B10/10wesn ";

open(CSH,">plot_movie.csh");
print CSH "gmtset BASEMAP_TYPE plain ANOT_FONT_SIZE 9 HEADER_FONT_SIZE 15\n";
print CSH "makecpt -Cpolar -T0/255/2 -Z -V > grd.cpt\n";

print "Processing locations: $xy_file\n";

open FILE, "$xy_file" or die $!;
@lines = <FILE>;
close(FILE);
$nlines = @lines;
print "  number of lines: $nlines\n\n";


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
  print "file min/max: $minmax \n  min: $min\n  max: $max\n";

  # transforms displacements to range [0,255] for colors
  if( abs($min) > abs($max) ){$max = abs($min);}
  open FILE, ">$file.xyz" or die $!;
  for($i=0;$i<$nlines;$i++){
    # scale between 0,1
    $val = ($lines_f[$i] + $max)/(2.0*$max);
    # scale between 0, 255
    $val = $val * 255.0;

    $coord = $lines[$i];
    chomp($coord);
    print FILE "$coord $val \n";
  }
  close(FILE);

  # determines new min/max of scaled displacements
  $minmax2 = `minmax $file.xyz | awk '{print \$7}' | sed -e "s:<: :" | sed -e "s:>: :" `;
  chomp($minmax2);
  ($min2,$max2) = split("/",$minmax2);
  print "  min2/max2: $minmax2 \n\n";

  #print CSH "paste $xy_file $file > $file.xyz \n";

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
