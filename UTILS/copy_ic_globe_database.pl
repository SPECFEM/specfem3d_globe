#!/usr/bin/perl -w

# this script collects the required mesh database files to plot IC in Paraview

use POSIX;

if (@ARGV != 3 and @ARGV != 4) {die("copy_ic_globe_database.pl slice_file lsf_machine_file filename [jobid]\n");}
$sfile = $ARGV[0];
$machine = $ARGV[1];
$filename = $ARGV[2];
if (@ARGV == 4) {$jobid = $ARGV[3];} else {$jobid ="";}

open(FILE,"$sfile") or die("Error opening file $sfile\n");
(@slices) = <FILE>;
close(FILE);

($nline) = split(" ",`wc $machine`);
if ($nline != 1) 
   {die("Check if machine file is LSF machine file format: if not, modify perl script\n");}
# LSF machine
open(FILE,"$machine") or die("Error opening file $machine\n");
($machines) = <FILE>;
close(FILE);
(@mac) = split(" ",$machines); $sl = 0;
for ($i=0;$i<@mac;$i=$i+2) {
  $sl = $sl + $mac[$i+1];
}
print "Total number of slices = $sl\n";
for ($i=0;$i<@mac;$i=$i+2) {
  $num_slices = $mac[$i+1];
  for ($j=0;$j<$num_slices;$j++) {
    $sl = $sl  - 1;
    $slice_to_node{$sl} = $mac[$i];
  }
}

for($i=0;$i<@slices;$i++) {
  ($slice[$i]) = split(" ",$slices[$i]);
  $node[$i] = $slice_to_node{$slice[$i]};
  print "$slice[$i], $node[$i]\n";
}

@files = ("AVS_DXelements","AVS_DXpoints","array_dims","$filename","solver_data_1");
@exts = ("txt","txt","txt","bin","bin");


for($i=0;$i<@slices;$i++) {
  for($j=0;$j<@files;$j++){
  $string = sprintf("scp $node[$i]:/scratch/$ENV{USER}/DATABASES_MPI.$jobid/proc%06d_reg3_$files[$j].$exts[$j] .", $slice[$i]);
  print "$string\n";
  system("$string");
}}

