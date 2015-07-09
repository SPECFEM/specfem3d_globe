#!/usr/bin/perl

# Mesh an inflated cube and write an opendx file for line rendering.
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
my $nbnode=20;
my $nbline=24;
my @iaddx = (0,1,2,2,2,1,0,0,0,1,2,2,2,1,0,0,0,2,2,0);
my @iaddy = (0,0,0,1,2,2,2,1,0,0,0,1,2,2,2,1,0,0,2,2);
my @iaddz = (0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,1,1,1,1);
my $radius = 100;
my @order=([0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,0],[8,9],[9,10],[10,11],[11,12],[12,13],[13,14],[14,15],[15,8],[0,16],[16,8],[2,17],[17,10],[4,18],[18,12],[6,19],[19,14]);
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
my $nblignes = ($nbpoints/20)*24;
open(OUT, ">mydxfile3D.dx");
print OUT " object 1 class array type float rank 1 shape 3 items         ",$nbpoints,"  data follows\n";
foreach (@points)
{ print OUT $_,"\n";
}
print OUT " object 2 class array type int rank 1 shape 2 items         ",$nblignes,"  data follows\n";
for (my $i=0;$i<$nb_elem_x*$nb_elem_y*$nb_elem_z;$i++)
{ for (my $j=0;$j<$nbline;$j++)
  { print OUT ($i*$nbnode+$order[$j]->[0]," ",$i*$nbnode+$order[$j]->[1],"\n");
  }
}
print OUT <<"EOF"
  attribute \"element type\" string \"lines\"
 attribute \"ref\" string \"positions\"
 object 3 class array type float rank 0 items $nblignes  data follows
EOF
;
for (my $k=0;$k<$nblignes;$k++)
{ print OUT (int($k/24),"\n");
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
