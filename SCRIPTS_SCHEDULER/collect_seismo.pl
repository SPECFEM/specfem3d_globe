#!/usr/bin/perl -w

if (@ARGV != 1) {die("collect_seismo.pl machinefile\n");}

$machine = $ARGV[0];

open(FILE,"$machine") or die("Error opening file $machine\n");
(@junk) = <FILE>;
close(FILE);

for($i=0;$i<@junk;$i++) {
  ($node) = split(" ",$junk[$i]);
  push(@nodes,$node);
}

#print "@nodes\n";
foreach $node (@nodes) {
    system("scp $node:/scratch/$ENV{USER}/DATABASES_MPI/*sem* .");
    print "$node\n";}

`shmux -M50 -Sall -c "rm -f /scratch/$ENV{USER}/DATABASES_MPI/*sem?" - < $machine`;
