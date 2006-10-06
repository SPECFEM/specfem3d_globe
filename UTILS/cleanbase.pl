#!/usr/bin/perl -w

if (@ARGV != 1) {die("cleanbase.pl machinefile\n");}
$mymachine = $ARGV[0];
if (not -f $mymachine) {die("check if $mymachine is a file or not\n");}

#`cluster-fork --pe-hostfile $mymachine --verbose "rm -rf /state/partition1/scratch/$ENV{USER}/DATABASES_MPI/*" `;
`shmux -M50 -Sall -c "rm -rf /state/partition1/scratch/$ENV{USER}/*; mkdir /state/partition1/scratch/$ENV{USER}/DATABASES_MPI;" - < $mymachine `;

