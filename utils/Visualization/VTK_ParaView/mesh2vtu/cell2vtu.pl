#!/usr/bin/perl

use Getopt::Std;
use FindBin;

$cell2vtu = "$FindBin::Bin/cell2vtu";

sub Usage {
    print STDERR <<END;
Usage: $0 -i input-file -o output-file
    Takes an input file (binary) with a number of points and a number of cells
    and transforms them into an unstructured grid file

    -i input-file (Binary file)
    -o output-file (XML Unstructured Grid File)

    Input Binary files have this structure:
      number_of_points          integer (4 bytes)
      x_1, y_1, z_1         3 reals (4 bytes each)
      ...
      x_n, y_n, z_n         3 reals (4 bytes each)
      number_of_cells           integer (4 bytes)
      cell_1 (eight points), scalar_1    8 integers (4 bytes each), 1 real (4 bytes)
      ...
      cell_n, scalar_n                   8 integers (4 bytes each), 1 real (4 bytes)

    This is a wrapper around cell2vtu

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
#print "$cell2vtu $opt_i $opt_o\n";
system("$cell2vtu $opt_i $opt_o");

1;
