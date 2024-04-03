#!/usr/bin/perl -w

# this script cleans ONLY the scratch directory given by Par_file
#   assumes that machine file lists the nodes line by line
#   requires 'shmux'
# Qinya Liu, Caltech, May 2007

if (@ARGV != 2) {die("cleanbase_lsf_multi.pl machinefile Par_file\n");}

$machine = $ARGV[0];
$par_file = $ARGV[1];

open(FILE3,"<$par_file") or die ("Fatal Error openning file $par_file\n");
while (<FILE3>) {
   if ($_ =~ /^LOCAL_PATH/) {
  chop;
  @vals = split("=", $_);
  $mpidir = $vals[1];
  $mpidir =~ s/^\s+//;
  $mpidir =~ s/\s+$//;
  close(FILE3);
  last;
   }
}

`shmux -M50 -Sall -c "rm -rf $mpidir" - < $machine`;

