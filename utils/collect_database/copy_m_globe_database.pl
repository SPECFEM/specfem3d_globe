#!/usr/bin/perl -w

# this script collects the required mantle database files to prepare mesh
# files for Paraview
#   array_dims.txt, solver_data_2.bin and kernel_filename/solver_data_1.bin
#   you need to supply the actual lsf_machine_file
#   (1 line, node name and number of procs used)
#   so that the correct machine file can be selected for the given slices (in slice_file)
# if you are collecting anything other than kernel files, solver_data_1.bin is
#   collected instead

# Qinya Liu, Caltech, May 2007


use POSIX;

if (@ARGV != 3 and @ARGV != 4) {die("copy_m_globe_database.pl slice_file lsf_machine_file filename [jobid]\n");}
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

@files = ("array_dims","solver_data_2","$filename");
@exts = ("txt","bin","bin");


for($i=0;$i<@slices;$i++) {
  for($j=0;$j<@files;$j++){
  $string = sprintf("scp $node[$i]:/scratch/$ENV{USER}/DATABASES_MPI$jobid/proc%06d_reg1_$files[$j].$exts[$j] .", $slice[$i]);
  print "$string\n";
  system("$string");
}}

