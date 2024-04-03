#!/usr/bin/perl

# Mesh an inflated cube and write an opendx file for hexahedron rendering.
# David Michea. University of Pau.
# March 2007.

print "\nnb d'elements radiaux ? : ";
$nbelem=<STDIN>;
chomp($nbelem);
print "\ndeformation (0->1) ? : ";
$deformation=<STDIN>;
chomp($deformation);


my $pi=3.141592653589793;
my $pisur2 = (3.141592653589793 / 2);
my $nb_elem_x=$nbelem;
my $nb_elem_y=$nbelem;
my $nb_elem_z=$nbelem;
my $nbnode=8;
my @iaddx = (0,0,2,2,0,0,2,2);
my @iaddy = (2,0,2,0,2,0,2,0);
my @iaddz = (2,2,2,2,0,0,0,0);
my $radius = 100;
my @points;
# calcul des coordonnees
for (my $ix=0;$ix<2*$nb_elem_x;$ix+=2)
{ for (my $iy=0;$iy<2*$nb_elem_y;$iy+=2)
  { for (my $iz=0;$iz<2*$nb_elem_z;$iz+=2)
    { for (my $inode=0;$inode<$nbnode;$inode++)
      { my ($x,$y,$z) = compute_coord($ix+$iaddx[$inode], $iy+$iaddy[$inode], $iz+$iaddz[$inode],$nb_elem_x*2,$nb_elem_y*2,$nb_elem_z*2,$radius);
        push (@points, join(" ",($x,$y,$z)));
      }
    }
  }
}
# ecriture des resultats dans fichier opendx
my $nbpoints = scalar(@points);
my $nbcubes = $nbpoints/$nbnode;
open(OUT, ">mydxfile3D.dx");
print OUT " object 1 class array type float rank 1 shape 3 items         ",$nbpoints,"  data follows\n";
foreach (@points)
{ print OUT $_,"\n";
}
print OUT " object 2 class array type int rank 1 shape 8 items         ",$nbcubes,"  data follows\n";
for (my $i=0;$i<$nb_elem_x*$nb_elem_y*$nb_elem_z;$i++)
{ for (my $j=0;$j<$nbnode;$j++)
  { print OUT " ",$i*$nbnode+$j;
  }
  print OUT "\n";
}
print OUT <<"EOF"
  attribute \"element type\" string \"cubes\"
 attribute \"ref\" string \"positions\"
 object 3 class array type float rank 0 items $nbcubes  data follows
EOF
;
for (my $k=1;$k<=$nbcubes;$k++)
{ print OUT 1,"\n";
}
print OUT <<EOF
 attribute "dep" string "connections"
 object "irregular positions irregular connections" class field
 component "positions" value 1
 component "connections" value 2
 component "data" value 3
EOF
;
close(OUT);

sub compute_coord
{ my ($ix,$iy,$iz,$nb_elem_x,$nb_elem_y,$nb_elem_z,$radius) = @_;
  my ($x,$y,$z);

  my $ratio_x = $ix/$nb_elem_x;
  my $ratio_y = $iy/$nb_elem_y;
  my $ratio_z = $iz/$nb_elem_z;

  my $factx = (2*$ratio_x -1);
  my $facty = (2*$ratio_y -1);
  my $factz = (2*$ratio_z -1);

  my $xi = $pisur2*$factx;
  my $eta = $pisur2*$facty;
  my $gamma = $pisur2*$factz;

  $x=($radius / sqrt(3)) * $factx  * (1 + (cos($eta)+cos($gamma))*$deformation / $pi);
        $y=($radius / sqrt(3)) * $facty  * (1 + (cos($xi)+cos($gamma))*$deformation / $pi);
  $z=($radius / sqrt(3)) * $factz  * (1 + (cos($xi)+cos($eta))*$deformation / $pi);

  return ($x,$y,$z);
}
