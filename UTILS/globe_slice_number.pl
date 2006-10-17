#!/usr/bin/perl -w
#  This program figures out the global slice number for given simulation files

# modify the following line for the correct location of perl libs (UTILS/lib)
use lib '/opt/seismo-util/lib/perl';
use CMT_TOOLS;

if (@ARGV != 3) {die("Usage: globe_slice_number.pl CMTSOLUTION STATIONS_ADJOINT Par_file\n");}
$cmt = $ARGV[0];
$sta = $ARGV[1];
$par_file = $ARGV[2];

$matlab_exe="matlab_r14";
$m_dir="matlab";
# obtain event location
($elat,$elon) = get_cmt_location($cmt);

# obtain station location
($slat,$slon) = split(" ",`awk 'NR == 2 {print \$3, \$4}' $sta`);

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

# check matlab program
$matlab = `which matlab >& /dev/null ; echo $?`; chomp($matlab);
if ($matlab) {die("Check if MATLAB program is installed in your system\n");}

# minor 
print "running matlab program ...\n";
open(MATLAB, "| $matlab_exe -nosplash -nojvm > matlab.txt");
print MATLAB "cd $m_dir \n";
if ($nchunk == 6) {
  print  MATLAB "slice_number($elat,$elon,$slat,$slon,$nproc,0)\n";
  print  "slice_number($elat,$elon,$slat,$slon,$nproc,0)\n";
}else {
  print  MATLAB "slice_number($elat,$elon,$slat,$slon,$nproc,0,$nchunk,$xi_width,$eta_width,$clon,$clat,$grot)\n";
 print  "slice_number($elat,$elon,$slat,$slon,$nproc,0,$nchunk,$xi_width,$eta_width,$clon,$clat,$grot)\n";
}

print MATLAB "quit\n";
close(MATLAB);
print "output to slices_minor file\n";
system(" awk 'NR >= 14 && NF > 0 {print \$0}' matlab.txt | grep -v '>' > slices_minor ");

if ($nchunk == 6) {
# major
open(MATLAB, "| $matlab_exe -nosplash -nojvm > matlab.txt");
print MATLAB "cd $m_dir \n";
print MATLAB "slice_number($elat,$elon,$slat,$slon,$nproc,1)\n";
print MATLAB "quit\n";
close(MATLAB);
print "output to slices_major file\n";
system(" awk 'NR >= 14 && NF > 0 {print \$0}' matlab.txt | grep -v '>' > slices_major ");

print "output to slices_ic file\n";
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
}

# figure out the normal to the source-receiver plane
print "calculate normal to source and receiver plane\n";
open(MATLAB, "| $matlab_exe -nosplash -nojvm > matlab.txt");
print MATLAB "cd $m_dir \n";
print  MATLAB "[s(1),s(2),s(3)] = normal_plane($elat,$elon,$slat,$slon);\ns\n";
print MATLAB "quit\n";
close(MATLAB);
system(" awk 'NR >= 14 && NF > 0 {print \$0}' matlab.txt | grep -v '>' > normal_plane.txt ");
$normal_plane = `awk 'NR >= 14 && NF > 0 {print \$0}' matlab.txt | grep -v '>' `;
print "$normal_plane \n";

# clean up temp files
system(" rm -f matlab.txt slice_ab_old");

