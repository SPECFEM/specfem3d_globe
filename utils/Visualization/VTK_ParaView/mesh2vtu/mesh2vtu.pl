#!/usr/bin/perl

use Getopt::Std;
use FindBin;

$mesh2vtu = "$FindBin::Bin/mesh2vtu";

sub Usage {
    print STDERR <<END;
Usage: $0 -i input-file -o output-file
    Takes an input file (binary) with a number of points and a number of cells
    and transforms them into an unstructured grid file

    -i input-file (Binary file)
    -o output-file (XML Unstructured Grid File)

    Input Binary files have this structure:
      number_of_points          integer (4 bytes)
      x_1, y_1, z_1, scalar_1   4 reals (4 bytes each)
      ...
      x_n, y_n, z_n, scalar_n   4 reals (4 bytes each)
      number_of_cells           integer (4 bytes)
      cell_1 (eight points)     8 integers (4 bytes each)
      ...
      cell_n                    8 integers (4 bytes each)

    This is a wrapper around mesh2vtu

    Brian Savage 6/26/2004
    
END
    exit(-1);
}

if(@ARGV == 0) {
    Usage ();
}

if(!getopts('i:o:')){die "Check input paramaters \n";}

if(!defined($opt_i)) {
    die "$0: Must specify input file -i input-file\n";
}
if(!defined($opt_o)) {
    die "$0: Must specify output file -o output-file\n";
}
#print "$mesh2vtu $opt_i $opt_o\n";
system("$mesh2vtu $opt_i $opt_o");

1;
