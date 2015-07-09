#!/usr/bin/perl

#
#  Script to change the version number in f90 codes
#
#  Author : Dimitri Komatitsch, EPS - Harvard University, May 1998
#

#
# read all f90 and F90 (and *.h) files in the current directory
# f90 files are supposed to have the extension "*.f90" or "*.F90" or "*.h"
#

#
# known bug : does the changes also in constant strings, and not only
# in the code (which is dangerous, but really easier to program...)
#

@objects = `ls src/*/*.f90 src/*/*.F90 src/*/*.c src/*/*.cu src/*/*.h.in src/*/*.h setup/*.h.in`;

foreach $name (@objects) {
  chop $name;

# change tabs to white spaces
  system("expand -2 < $name > _____tutu01_____");
  $f90name = $name;
  print STDOUT "Changing word in file $name ...\n";

  open(FILEF77,"<_____tutu01_____");
  open(FILEF90,">$f90name");

# open the source file
  while($line = <FILEF77>) {
    chop $line;

# suppress trailing white spaces and carriage return
    $line =~ s/\s*$//;

# converts tabs to space
#      $line =~ s/\t/ /g;

# change the version number
#    $line =~ s#!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1#!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0#og;

# change the copyright information

#      $line =~ s#!             and University of Pau / CNRS / INRIA, France#!             and CNRS / INRIA / University of Pau, France#og;
#      $line =~ s#! \(c\) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA#! \(c\) Princeton University and CNRS / INRIA / University of Pau#og;
#     $line =~ s#rmass_sigma#rmass_time_integral_of_sigma#og;

# write the modified line to the output file
    print FILEF90 "$line\n";
  }

  close(FILEF77);
  close(FILEF90);
}

system("rm -f _____tutu01_____");

