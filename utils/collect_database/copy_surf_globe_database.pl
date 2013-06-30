#!/usr/bin/perl -w

# this script collects the required database (kernel file or solver_data_1.bin)
# to plot surface mesh in Paraview
# Qinya Liu, May 2007, Caltech

use POSIX;

if (@ARGV != 3 and @ARGV != 4) {die("copy_surf_globe_database.pl slice_file lsf_machine_file filename [jobid]\n");}
$sfile = $ARGV[0];
$machine = $ARGV[1];
$filename = $ARGV[2];

if ($filename !~ /kernel/) {$filename="solver_data_1";}
if (@ARGV == 4) {$jobid = ".$ARGV[3]";} else {$jobid ="";}

open(FILE,"$sfile") or die("Error opening file $sfile\n");
(@slices) = <FILE>;
close(FILE);
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

@regs=("reg1","reg1","reg1","reg1","reg1","reg2","reg2","reg2","reg2");
@files = ("array_dims","$filename","boundary_disc","boundary","solver_data_2",  "array_dims","$filename","boundary","solver_data_2");
@exts = ("txt","bin","bin","bin","bin","txt","bin","bin","bin");


for($i=0;$i<@slices;$i++) {
  for($j=0;$j<@files;$j++){
  $string = sprintf("scp $node[$i]:/scratch/$ENV{USER}/DATABASES_MPI$jobid/proc%06d_$regs[$j]_$files[$j].$exts[$j] .", $slice[$i]);
  print "$string\n";
  system("$string");
}}

