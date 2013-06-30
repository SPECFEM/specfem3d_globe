#!/usr/bin/perl -w

# completely cleans your scratch disk, and regenerate DATABASES_MPI/ directory in the scratch disk
#   uses 'shmux' to have simultaneous access to all nodes

# Qinya Liu, Caltech, May 2007


if (@ARGV != 1) {die("cleanbase.pl machinefile\n");}
$mymachine = $ARGV[0];
if (not -f $mymachine) {die("check if $mymachine is a file or not\n");}

`shmux -M50 -Sall -c "rm -rf /scratch/$ENV{USER}/*; mkdir /scratch/$ENV{USER}/DATABASES_MPI;" - < $mymachine `;

