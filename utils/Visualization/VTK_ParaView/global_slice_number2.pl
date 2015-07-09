#!/usr/bin/perl -w

#  This program figures out the global slice number on a source-receiver
#  belt for given simulation parameters
#  uses xglobal_slice_number and xnormal_plane, compile first from global_slice_util/
#  Qinya Liu, Caltech, May 2007

# modify the following line for the correct location of perl libs (utils/lib)
use FindBin;
use lib "$FindBin::Bin/../../lib";
use CMT_TOOLS;

if (@ARGV != 4) {die("Usage: global_slice_number.pl CMTSOLUTION STATIONS_ADJOINT Par_file lat_belt_width (in degrees)\n");}
$cmt = $ARGV[0];
$sta = $ARGV[1];
$par_file = $ARGV[2];
$lat0 = $ARGV[3];

# obtain event location
($elat,$elon) = get_cmt_location($cmt);

# obtain station location
($slat,$slon) = split(" ",`awk 'NR == 1 {print \$3, \$4}' $sta`);

print "event = ($elat,$elon);  station = ($slat,$slon) \n";

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

# minor
print "compute slices along minor arc ...\n";
if ($nchunk == 6) {
  $command = "xglobal_slice_number $elon $elat $slon $slat $nproc 0 $lat0";
}else {
  $command = "xglobal_slice_number $elon $elat $slon $slat $nproc 0 $nchunk $xi_width $eta_width $clon $clat $grot $lat0";}
system(" $command > slices_minor ");

# major
print "compute slices along major arc ...\n";
if ($nchunk == 6) {
  $command = "xglobal_slice_number $elon $elat $slon $slat $nproc 1 $lat0";
}else {
  $command = "xglobal_slice_number $elon $elat $slon $slat $nproc 1 $nchunk $xi_width $eta_width $clon $clat $grot $lat0";}
system(" $command > slices_major ");

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
$command = "xnormal_plane $elat $elon $slat $slon";
$result = `$command`;
print "$result\n";

# clean up temp files
system(" rm -f slice_ab_old gcarc_station.txt xsection_translate.txt");

