#!/usr/bin/perl
#
# reads in binary files produced from xcreate_movie_GMT_global
# and creates png files with transparency
#
# uses segmented M-file to plot with psxy
#

use Getopt::Std;
use POSIX;

sub Usage {
  print STDERR <<EOF;

Usage: plot_movie_GMT_binary_PNG.pl file_name

  ex. ./plot_movie_GMT_binary_PNG.pl OUTPUT_FILES/bin_movie_*.d

EOF
exit(1)
}

if (!getopts('l:L:cdts')) {die('Check input arguments\n');}
@ARGV > 0 or Usage();


#######################################################
## PARAMETERS

# point locations
$xy_file = "OUTPUT_FILES/bin_movie.xy";

# specfem3D Par_file settings
$nex = 160;
$nproc = 5;
$nchunks = 6;

# non-linear scaling
$power_scaling = 0.5;

# global region
$R = "-Rd";

#plate carre projection
$JM = "-JQ0/0/15";

# color tables
$colorcpt = "blue_white_red.cpt";
$maskcpt = "gray_pyramid_inv.cpt";

#######################################################

# verbosity
$verbose = ""; # "-V";
# grid
$grdfile = "movie.grd";
# annotation
$B = " -B10/10wesn ";


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


open(CSH,">plot_movie.csh");
print CSH "gmtset BASEMAP_TYPE plain ANOT_FONT_SIZE 9 HEADER_FONT_SIZE 15\n";

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

  # restore gzipped file
  if( $gzipped == 1 ) {system("mv $filegz.org $filegz");}

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

    # calculates a nice rectangle size
    ($lon,$lat) = split(" ",$coord);
    $rect_length = abs($lat) / 90. ;
    $rect_length = 0.03 + ($rect_length**2) * 0.1;
    $rect_height = 0.03 + ($rect_length**2) * 0.005;

    # adds point
    print FILE "$coord $val $rect_length $rect_height \n";

    # stores values in arrays
    push(@elem_val, $val );
    $line = "$coord $val ";
    push(@elem_corners, $line);
  }
  close(FILE);


  # specfem3D output for MOVIE_COARSE
  $nex_per_proc = $nex / $nproc  ;
  $total = $nproc * $nex_per_proc * $nproc * $nex_per_proc * $nchunks;
  if( $total != $nlines ){ die("error nex/proc $nlines $total");}

  # writes out segment file
  open MFILE, ">$file.M.xyz" or die $!;
  for($n=0;$n<$nproc*$nproc*$nchunks;$n++){
    for($m=0;$m<$nex_per_proc-1;$m++){
      for ($k=0;$k<$nex_per_proc-1;$k++){
        # global array index
        $i = $k + $m*$nex_per_proc + $n*$nex_per_proc*$nex_per_proc;
        $val = $elem_val[$i];
        # adds element
        print MFILE "> -Z$val -W- \n";
        print MFILE "$elem_corners[$i] \n";
        print MFILE "$elem_corners[$i+1] \n";
        print MFILE "$elem_corners[$i+$nex_per_proc+1] \n";
        print MFILE "$elem_corners[$i+$nex_per_proc] \n";
      }
    }
  }
  close(MFILE);

  $ps_file = "$file.ps";

  # start ps-file
  print CSH "psxy $file.xyz $JM $R $B -Sr -C$colorcpt -K -P $verbose > $ps_file\n";
  print CSH "psxy $file.M.xyz $JM $R $B -M -L -N -C$colorcpt -P -O  $verbose >> $ps_file\n";
  # end ps-file

  # start ps-file mask
  print CSH "psbasemap $JM $R $B -G0 -P $verbose -K > $ps_file.mask\n";
  print CSH "psxy $file.xyz $JM $R $B -Sr -C$maskcpt -P -K -O $verbose >> $ps_file.mask\n";
  print CSH "psxy $file.M.xyz $JM $R $B -M -L -N -C$maskcpt -P -O $verbose >> $ps_file.mask\n";
  # end ps-file mask

  # convert to png
  print CSH "convert -crop 1080x540+6+6 $ps_file $file.disp.png\n";
  print CSH "convert -crop 1080x540+6+6 $ps_file.mask $file.mask.png\n";

  # creates file with transparency (opacity)
  print CSH "composite -compose CopyOpacity $file.mask.png $file.disp.png $file.png\n";

  # remove temporary data file
  print CSH "rm -f $file.xyz $file.M.xyz\n";

  # removes temporary images
  print CSH "rm -f $ps_file $ps_file.mask  \n";
  print CSH "rm -f $file.mask.png $file.disp.png \n";

  print CSH "echo \n";
  print CSH "echo 'plotted: $file.png' \n";
  print CSH "echo \n";

}
close(CSH);

print "\nplotting... \n\n";

system("csh -f plot_movie.csh");

system("rm -f plot_movie.csh");
system("rm -f $grdfile $grdfile.1");


