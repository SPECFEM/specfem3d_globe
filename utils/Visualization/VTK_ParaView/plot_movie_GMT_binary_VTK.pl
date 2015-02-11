#!/usr/bin/perl -W
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

Usage: $0 file_name

  ex. $0 OUTPUT_FILES/bin_movie_00***.d

EOF
exit(1)
}

if (!getopts('l:L:cdts')) {die('Check input arguments\n');}
@ARGV > 0 or Usage();


#######################################################
## PARAMETERS

# specfem3D Par_file settings
$nex = 160;
$nproc = 5;
$nchunks = 6;

# gunzipped files (ending .gz)
$gzipped = 0;

# point locations
$xy_file = "OUTPUT_FILES/bin_movie.xy";

# non-linear scaling
$power_scaling = 0.5;


#######################################################


# gunzips location file
if( $gzipped == 1 ){
  system("cp -v $xy_file.gz $xy_file.gz.org");
  system("gunzip -f $xy_file.gz");
}

print "Processing locations: $xy_file\n";

(@lines) = read_binary_file_locations($xy_file);
$nlines = @lines;
print "  number of lines: $nlines\n\n";

# restore gzipped file
if( $gzipped == 1 ) {system("mv $xy_file.gz.org $xy_file.gz");}

open(CSH,">plot_movie.csh");


foreach $file (@ARGV) {

  if (not -f $file) {die("No $file\n");}

  print "Processing frame $file...\n";

  # gunzips data file
  if( $gzipped == 1 ){
    $filegz = $file;
    $file = substr($file,0,length($filegz)-3);
    system("cp $filegz $filegz.org");
    system("gunzip -f $filegz");
  }
  print "    $file\n";

  # reads displacement file (binary data)
  (@lines_f) = read_binary_file_data($file);
  $nlines_f = @lines_f;

  # restore gzipped file
  if( $gzipped == 1 ) {system("mv $filegz.org $filegz");}

  if ($nlines_f != $nlines) {die("number of lines differ\n");}


  # non-linear scaling of displacement data and min/max
  $min = 1.e30;
  $max = -1.e30;
  for($i=0;$i<$nlines;$i++){
    # scale between 0,1
    $val = $lines_f[$i];
    if( $power_scaling > 0 ) {
      if( $val > 0 ){
        $val = $val**$power_scaling;
      }else{
        $val = - abs($val)**$power_scaling;
      }
    }
    # determines min/max of displacements
    if( $val < $min) {$min = $val;}
    if( $val > $max) {$max = $val;}
    $lines_f[$i] = $val;
  }
  if( abs($min) > abs($max) ){$max = abs($min);}

  # writes out segment file
  open VTKFILE, ">$file.vtk" or die $!;
  print VTKFILE "# vtk DataFile Version 3.1\n";
  print VTKFILE "specfem3D_data\n";
  print VTKFILE "ASCII\n";
  print VTKFILE "DATASET POLYDATA\n";
  print VTKFILE "POINTS $nlines float\n";

  @elem_val = ();
  for($i=0;$i<$nlines;$i++){
    # scale between 0,1
    $val = ($lines_f[$i] + $max)/(2.0*$max);
    # scale between 0, 255
    $val = $val * 255.0;

    # prints value together with coordinates
    $coord = $lines[$i];
    chomp($coord);

    # adds point: uses lon, lat and 0 (flat earth file)
    print VTKFILE "$coord 0.0 \n";

    # stores values in arrays
    push(@elem_val, $val );
  }
  print VTKFILE "\n";


  # specfem3D output for MOVIE_COARSE
  $nex_per_proc = $nex / $nproc  ;
  $total = $nproc * $nex_per_proc * $nproc * $nex_per_proc * $nchunks;
  if( $total != $nlines ){ die("error nex/proc $nlines $total");}

  $npoly = $nproc * $nproc* $nchunks * ($nex_per_proc-1) * ($nex_per_proc-1);
  print VTKFILE "POLYGONS $npoly ",$npoly*5," \n";
  $count = 0;
  for($n=0;$n<$nproc*$nproc*$nchunks;$n++){
    for($m=0;$m<$nex_per_proc-1;$m++){
      for ($k=0;$k<$nex_per_proc-1;$k++){
        # global array index
        $i = $k + $m*$nex_per_proc + $n*$nex_per_proc*$nex_per_proc;

        # adds element
        print VTKFILE "4 ",$i," ",$i+1," ",$i+$nex_per_proc+1," ",$i+$nex_per_proc," \n";
        $count++;
      }
    }
  }
  print VTKFILE "\n";

  print "count: $count \n";

  print VTKFILE "POINT_DATA $nlines \n";
  print VTKFILE "SCALARS displacement float\n";
  print VTKFILE "LOOKUP_TABLE default\n";
  for($i=0;$i<$nlines;$i++){
    $val = $elem_val[$i];
    print VTKFILE "$val \n";
  }
  print VTKFILE "\n";
  close(VTKFILE);

  # creates png files
  #print CSH "python plot_VTKdisp.py $file.vtk \n";
  #print CSH "python plot_VTKdisp_gray.py $file.vtk \n";
  #print CSH "mv bin_disp.png $file.disp.png \n";
  #print CSH "mv bin_mask.png $file.mask.png \n";

  # creates file with transparency (opacity)
  #print CSH "composite -compose CopyOpacity $file.mask.png $file.disp.png $file.png\n";

  print CSH "echo \n";
  print CSH "echo 'plotted: $file.png' \n";
  print CSH "echo \n";

}
close(CSH);

print "\nplotting... \n\n";

system("csh -f plot_movie.csh");

system("rm -f plot_movie.csh");

#-------------------------------------------------------------

sub read_binary_file_locations{

  my($xy_file) = @_;
  my(@return_lines);
  my(@lines,$line,$is_ok,$junk,$n1,$n2,$data1,$data2);

  # reads in locations (binary data)
  open FILE, "$xy_file" or die $!;
  binmode FILE;

  @lines= () ;
  $line="";
  $is_ok = 1;
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

  }
  close(FILE);

  @return_lines = @lines;
  return @return_lines;

}


#-------------------------------------------------------------

sub read_binary_file_data{

  my($file) = @_;
  my(@return_lines_f);
  my(@lines_f,$is_ok,$junk,$n1,$data1);

  open FILE, "$file" or die $!;
  binmode FILE;

  @lines_f = () ;
  $line = "";
  $is_ok = 1;
  while ( $is_ok != 0) {
    read FILE, $junk, 4;

    $n1 = read FILE, $data1, 4;
    if( $n1 != 0 ){
      ($val) = unpack( "f",$data1);
    }

    read FILE, $junk, 4;

    #print "$n1 bytes read: value=$val \n";

    if( $n1 != 0 ){
      $line = "$val \n";
      push(@lines_f, $line);
      $is_ok = 1;
    }
    else{ $is_ok = 0;}
  }
  close(FILE);

  @return_lines_f = @lines_f;
  return @return_lines_f;

}

#-------------------------------------------------------------
