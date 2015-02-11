#!/usr/bin/perl -w

# this script collects the seismograms from the scratch directory given by Par_file on
# the machines given by machinefile

# Qinya Liu, Caltech, May 2007

if (@ARGV != 2) {die("collect_seismo_lsf_multi.pl machinefile Par_file\n");}

$machine = $ARGV[0];
$par_file = $ARGV[1];

# get the machine list
open(FILE,"$machine") or die("Error opening file $machine\n");
(@junk) = <FILE>;
close(FILE);

for($i=0;$i<@junk;$i++) {
  ($node) = split(" ",$junk[$i]);
  push(@nodes,$node);
}

# now get the LOCAL_PATH
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

foreach $node (@nodes) {
    system("scp $node:$mpidir/*sem* .");
    print "$node\n";}

# you can choose to delete the seismograms on the scratch disk after collecting them
#`shmux -M50 -Sall -c "rm -f $mpidir/*sem*" - < $machine`;

