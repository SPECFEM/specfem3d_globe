#! /usr/bin/perl

# count NGLL points in an OpenDX mesh file

# David Michea, University of Pau, France, January 2007

use Data::Dumper;

if (scalar(@ARGV) < 1)
{     print "usage : $0 mesh2analyse.dx [NGLL]\n";
      exit;
}

$filein = $ARGV[0];
$NGLL = $ARGV[1] || 0;
exit("file $filein does not exists") unless (-f $filein);

$numnode=0;
$nb_corners=8;
# opendx numbering of cubes' faces
my @faces_opendx = ([2,3,7,6], [0,2,6,4], [0,1,5,4], [1,3,7,5], [0,1,3,2], [4,5,7,6]);

# read datas
open (DF, "$filein");
while(my $line=<DF>)
{     chomp $line;
      if ($line =~ m/^.*object\s+1.*\s+(\d+)\s+data\s+follows\s*$/)
      {     $nbnode = $1;
      }
      if ($line =~ m/^(\s+[-]?\d+\.?\d*\s*)+$/)
      {     if ($line =~ m/\./) # use -?[0-9]+.[0-9]* notation for point's coordinate (ie : don't use scientific notation, not implemented yet)
            {     my @temp=split(/\s+/, $line);
                  shift @temp;
                  $nodes[$numnode] = \@temp;
                  $numnode++;
            }
            elsif ($line =~ m/(\s+\d+\s*){8}/)
            {     my @temp=split(/\s+/,$line);
                  shift @temp;
                  $hex[$numhex] = \@temp;
                  $numhex++;
            }
            else
            {     # couleur
            }
      }
}
close (DF);
exit "data read error  : use -?[0-9]+.[0-9]* notation for point's coordinate (ie : don't use scientific notation, not implemented yet)\n" unless ($nbnode==$numnode);
# dump datas into structure
my @cubes;
foreach $cube (@hex)
{     my @corners;
      foreach $corner (@{$cube})
      {     my @coords = @{$nodes[$corner]};
            push(@corners, \@coords);
      }
      push (@cubes, \@corners);
}
print "nb elements : $numhex\n";
my $faces = get_connexion_faces(\@cubes);
$nbfaces = scalar(@$faces);
print "nb de surfaces de contact entre ces elements : ",$nbfaces,"\n";
my $edges = get_connexion_edges($faces);
$faces=undef;
$edges = edge_clean_up($edges);
$nbedges = scalar(@$edges);
print "nb d'aretes de contact entre ces surfaces : ",$nbedges,"\n";
my $points = get_connexion_points($edges);
$edges=undef;
$points = point_clean_up($points);
$nbpoints = scalar(@$points);
$points=undef;
print "nb de points de contact entre ces aretes : ",$nbpoints,"\n";
print "\nNGLOB = $numhex.NGLL^3 - $nbfaces.NGLL^2 + $nbedges.NGLL - $nbpoints\n";
$res = (($numhex*($NGLL**3))-($nbfaces*($NGLL**2))+($nbedges*$NGLL)-$nbpoints);
print ("pour NGLL = $NGLL, NGLOB = $res\n") if ($NGLL);
exit;

sub point_clean_up
{     my $points_in = shift;
      my @points_out;
      while ($pt1=shift(@$points_in))
      {     my $suppr = 0;
            foreach $pt2 (@$points_in)
            {     if ($pt1->[0] == $pt2->[0] and $pt1->[1] == $pt2->[1] and $pt1->[2] == $pt2->[2])
                  {     $suppr++;
                        last;
                  }
            }
            push (@points_out, $pt1) unless ($suppr);
      }
      return \@points_out;
}

sub edge_clean_up
{     my $edges_in = shift;
      my @edges_out;
      while ($edge1=shift(@$edges_in))
      {     my $suppr = 0;
            foreach $edge2 (@$edges_in)
            {     if (eqedge($edge1, $edge2))
                  {     $suppr++;
                        last;
                  }
            }
            push (@edges_out, $edge1) unless ($suppr);
      }
      return \@edges_out;
}

sub get_connexion_points
{     my $edges = shift;
      my @conn_points;
      my $nb=0;
      foreach $couple (@$edges)
      {     $nb++;
            for (my $i=$nb;$i<$nbedges;$i++)
            {     my $tmp = get_conn_point($couple, $edges->[$i]);
                  push (@conn_points, $tmp) if ($tmp);
            }
      }
      return \@conn_points;
}

sub get_conn_point
{     my ($couple_1, $couple_2) = @_;
      foreach $pt1 (@$couple_1)
      {     foreach $pt2 (@$couple_2)
            {     return $pt1 if ($pt1->[0] == $pt2->[0] and $pt1->[1] == $pt2->[1] and $pt1->[2] == $pt2->[2]);
            }
      }
      return 0;
}

sub get_connexion_edges
{     my $faces = shift;
      my @conn_edges;
      my $nb=0;
      foreach $quad (@$faces)
      {     $nb++;
            my $cur_face_edges = get_face_edges($quad);
            for (my $i=$nb;$i<$nbfaces;$i++)
            {     my $tmp_face_edges = get_face_edges($faces->[$i]);
                  my $tmp = get_conn_edges($cur_face_edges, $tmp_face_edges);
                  push (@conn_edges, $tmp) if ($tmp);
            }
      }
      return \@conn_edges;
}

sub get_face_edges
{     my $quad = shift;
      my @edges;
      for(my $i=0;$i<4;$i++)
      {     my @edge_end;
            push (@edge_end, $quad->[$i-1]);
            push (@edge_end, $quad->[$i]);
            push (@edges, \@edge_end);
      }
      return \@edges;
}

sub get_conn_edges
{     my ($edges_set_1, $edges_set_2) = @_;
      foreach $edge1 (@$edges_set_1)
      {     foreach $edge2 (@$edges_set_2)
            {     if (eqedge($edge1,$edge2))
                  {     return $edge1;
                  }
            }
      }
      return 0;
}

sub eqedge
{     my ($edge1, $edge2) = @_;
      my $count = 0;
      foreach $pt1 (@$edge1)
      {     foreach $pt2 (@$edge2)
            {     $count++ if ($pt1->[0] == $pt2->[0] and $pt1->[1] == $pt2->[1] and $pt1->[2] == $pt2->[2]);
            }
      }
      return 1 if ($count == 2);
      return 0;
}

sub get_connexion_faces
{     my $cubes = shift;
      my @conn_faces;
      my $nb=0;
      foreach $hex (@$cubes)
      {     $nb++;
            my $cur_cube_faces = get_cube_faces($hex);
            for (my $i=$nb;$i<$numhex;$i++)
            {     my $tmp_cube_faces = get_cube_faces($cubes->[$i]);
                  my $tmp = get_conn_faces($cur_cube_faces, $tmp_cube_faces);
                  push (@conn_faces, $tmp) if ($tmp);
            }
      }
      return \@conn_faces;
}

sub get_cube_faces
{     my $hex = shift;
      my @faces;
      foreach my $face (@faces_opendx)
      {     my @face_corners;
            foreach (@$face)
            {     push (@face_corners, $hex->[$_]);
            }
            push (@faces, \@face_corners);
      }
      return \@faces;
}

sub get_conn_faces
{     my ($faces_set_1, $faces_set_2) = @_;
      foreach $face1 (@$faces_set_1)
      {     foreach $face2 (@$faces_set_2)
            {     if (eqface($face1,$face2))
                  {     return $face1;
                  }
            }
      }
      return 0;
}

sub eqface
{     my ($face1, $face2) = @_;
      my $count = 0;
      foreach $pt1 (@$face1)
      {     foreach $pt2 (@$face2)
            {     $count++ if ($pt1->[0] == $pt2->[0] and $pt1->[1] == $pt2->[1] and $pt1->[2] == $pt2->[2]);
            }
      }
      return 1 if ($count == 4);
      return 0;
}
