#!/usr/bin/perl

#
#  Script to change paths to home directory in the entire source code
#
#  Author : Dimitri Komatitsch, Caltech, May 2002, for Specfem3D Globe
#

#
# read all f90 files in the current directory
# f90 files are supposed to have the extension "*.f90"
#

      @objects = `ls *.f90`;

      foreach $name (@objects) {
            chop $name;

# copy input file
            system("cp $name __tempfile");
            $f90name = $name;
            print STDOUT "Changing paths to home directory in file $name ...\n";

            open(FILE_OLD_PATHS,"<__tempfile");
            open(FILE_NEW_PATHS,">$f90name");

# read f90 source file
      while($line = <FILE_OLD_PATHS>) {
            chop $line;

# suppress trailing white spaces and carriage return
      $line =~ s/\s*$//;

#
# change path to input data files in home directory
#
      $line =~ s#DATA#/S/home010/DATA#og;

#
# change path to output files produced by the code in home directory
# do not change path to include files used when compiling the code
#
      if(index($line,"include") < 0) { $line =~ s#OUTPUT_FILES#/S/home010/OUTPUT_FILES#og; }


# write this line
      print FILE_NEW_PATHS "$line\n";

      }

      close(FILE_OLD_PATHS);
      close(FILE_NEW_PATHS);

      }

# suppress temporary file
      system("rm -f __tempfile");

