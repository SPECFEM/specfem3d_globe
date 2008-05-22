#!/usr/bin/perl -w

if (@ARGV != 1) {die("remap_lsf_machines.pl machinefile\n");}

$machine = $ARGV[0];

open(FILE,"$machine") or die("Error opening file $machine\n");
(@junk) = <FILE>;
close(FILE);

for($i=0;$i<@junk;$i++) {
  @node_array = split(" ",$junk[$i]);
  foreach $node (@node_array) {
	next if ( $node =~ /^[0-9]/ );
  	push(@nodes, $node);
  }
}
foreach $node (@nodes) {
    print "$node\n";
}
