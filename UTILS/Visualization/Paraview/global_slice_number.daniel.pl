#!/usr/bin/perl -w

#  This program figures out the global slice number for given simulation parameters
#  uses xglobal_slice_number and xnormal_plane, compile first from global_slice_util/
# Qinya Liu, May 2007, Caltech

# modify the following line for the correct location of perl libs (UTILS/lib)
use lib '/opt/seismo-util/lib/perl';
use CMT_TOOLS;

if (@ARGV != 3) {die("Usage: global_slice_number.pl CMTSOLUTION STATIONS_ADJOINT Par_file\n");}
$cmt = $ARGV[0];
$sta = $ARGV[1];
$par_file = $ARGV[2];

# obtain event location
($elat,$elon) = get_cmt_location($cmt);

# obtain station location
@lat = ();
@lon = ();
open(IN,$sta);
@stations = <IN>;
close(IN);
$nstations = @stations;

for($i=0;$i<$nstations;$i++){
  $line = $stations[$i];
  (undef,undef,$slat,$slon,undef,undef) = split(" ",$line);
  push(@lat,$slat);
  push(@lon,$slon);

  print "event = ($elat,$elon);  station = ($lat[$i],$lon[$i]) \n";

}
#($slat,$slon) = split(" ",`awk 'NR == 1 {print \$3, \$4}' $sta`);
#print "event = ($elat,$elon);  station = ($slat,$slon) \n";

# obtain simulation parameters
if (not -f $par_file) {die("Check if $par_file exists or not\n");}
($nchunk) = split(" ",`grep NCHUNKS $par_file|awk '{print \$3}'`);
($nproc) = split(" ",`grep NPROC_XI $par_file|awk '{print \$3}'`);
print "NCHUNKS = $nchunk, nproc = $nproc\n";
if ($nchunk != 6) {
  ($xi_width) = split(" ",`grep ANGULAR_WIDTH_XI_IN_DEGREES $par_file|awk '{print \$3}'`);
  ($eta_width) = split(" ",`grep ANGULAR_WIDTH_ETA_IN_DEGREES $par_file|awk '{print \$3}'`);
  ($clat) = split(" ",`grep CENTER_LATITUDE_IN_DEGREES $par_file|awk '{print \$3}'`);
  ($clon) = split(" ",`grep CENTER_LONGITUDE_IN_DEGREES $par_file|awk '{print \$3}'`);
  ($grot) = split(" ",`grep GAMMA_ROTATION_AZIMUTH $par_file|awk '{print \$3}'`);
  print "xi_width = $xi_width, eta_width = $eta_width; clat = $clat, clon = $clon; grot = $grot\n";
}

# all
open(SLICE_ALL,">slices_all");
for ($i = 0; $i < $nproc*$nproc*$nchunk; $i++ ) {
    print SLICE_ALL "$i\n";
}
close(SLICE_ALL);  

# minor 
print "compute slices along minor arc ...\n";

for($istation=0;$istation<$nstations;$istation++){
  $slat = $lat[$istation];
  $slon = $lon[$istation];
  
  if ($nchunk == 6) {
    $command = "./xglobal_slice_number $elon $elat $slon $slat $nproc 0";
  }else {
    $command = "./xglobal_slice_number $elon $elat $slon $slat $nproc 0 $nchunk $xi_width $eta_width $clon $clat $grot";}
  
  if( $istation == 0 ){
    system(" $command > slices_minor ");
  }else {
    system(" $command >> slices_minor ");  
  }
}

#sorts out doubles
open(IN,"slices_minor");
@slices = <IN>;
close(IN);
$nslices = @slices;
@uni = ();
for($i=0;$i<$nslices;$i++){
  $sl = $slices[$i];
  $found = 0;
  for($j=0;$j<@uni;$j++){
    if( $sl == $uni[$j] ){ 
      $found = 1; 
      last;
    }
  }  
  if( $found == 0 ){
    push(@uni,$sl);
  }
}
open(OUT,">slices_minor");
for($i=0;$i<@uni;$i++){
  print OUT "$uni[$i]";
}
close(OUT);

# major
print "compute slices along major arc ...\n";

for($istation=0;$istation<$nstations;$istation++){
  $slat = $lat[$istation];
  $slon = $lon[$istation];
  
  if ($nchunk == 6) {
    $command = "./xglobal_slice_number $elon $elat $slon $slat $nproc 1";
  }else {
    $command = "./xglobal_slice_number $elon $elat $slon $slat $nproc 1 $nchunk $xi_width $eta_width $clon $clat $grot";}

  if( $istation == 0 ){
    system(" $command > slices_major ");
  }else {
    system(" $command >> slices_major ");  
  }
}

#sorts out doubles
open(IN,"slices_major");
@slices = <IN>;
close(IN);
$nslices = @slices;
@uni = ();
for($i=0;$i<$nslices;$i++){
  $sl = $slices[$i];
  $found = 0;
  for($j=0;$j<@uni;$j++){
    if( $sl == $uni[$j] ){ 
      $found = 1; 
      last;
    }
  }  
  if( $found == 0 ){
    push(@uni,$sl);
  }
}
open(OUT,">slices_major");
for($i=0;$i<@uni;$i++){
  print OUT "$uni[$i]";
}
close(OUT);





if ($nchunk == 6) {
  print "output slices_ic file\n";
  $nslice = $nproc * $nproc;
  `cat slices_minor slices_major | sort -g | awk '\$1 < $nslice {print \$1}' > slice_ab_old`;

  open(SLICE,"slice_ab_old");
  @tmp=<SLICE>;
  close(SLICE);

  open(SLICE_IC,">slices_ic");
  for ($i = 0; $i < $nslice; $i ++ ) {
    for ($j = 0; $j < @tmp; $j ++ ){
      if ($i == $tmp[$j]) {last;}
      if ($j == @tmp - 1) {print SLICE_IC "$i\n";}
    }
  }
  close(SLICE_IC);
}

# figure out the normal to the source-receiver plane
print "calculate normal to source and receiver plane\n";
$command = "./xnormal_plane $elat $elon $slat $slon";
$result = `$command`;
print "$result\n";

# clean up temp files
system(" rm -f slice_ab_old gcarc_station.txt xsection_translate.txt");

