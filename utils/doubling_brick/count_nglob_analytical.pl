#! /usr/bin/perl

# count NGLL points in an OpenDX mesh file (written initially for the doubling brick)

# -> limitations :  - brick have to be a parallelepiped having one edge parallel to (Ox)
#                   - points coordinates have to be in [-$hugeval; $hugeval]
#                   - and other ones ...

# be carefull : script one-shot written for a particular geometry, may have to be adapted, use it at your own risk ...

# David Michea, University of Pau, France, February 2007

if (scalar(@ARGV) < 1)
{     print "usage : $0 mesh2analyse.dx [NGLL]\n";
      exit;
}

$filein = $ARGV[0];
$NGLL = $ARGV[1] || 0;
exit("file $filein does not exists") unless (-f $filein);
$hugeval=1000000;
$numnode=0;
$nb_corners=8;
# opendx numbering of cube's faces
my @faces_opendx = ([2,3,7,6], [0,2,6,4], [0,1,5,4], [1,3,7,5], [0,1,3,2], [4,5,7,6]);

# read datas
open (DF, "$filein");
$phase=1;
while(my $line=<DF>)
{     chomp $line;
      if ($line =~ m/^.*object\s+1.*\s+(\d+)\s+data\s+follows\s*$/)
      {     $nbline = $1;
            if($phase!=1)
            {   print "error while reading data\n";
                exit;
            }
            for (my $i=0;$i<$nbline;$i++)
            {     my $line=<DF>;
                  chomp($line);
                  my @temp=split(/\s+/, $line);
                  shift @temp if (scalar(@temp)>3);
                  $nodes[$numnode] = \@temp;
                  $numnode++;
            }
            $phase=2;
      }
      elsif ($line =~ m/^.*object\s+2.*\s+(\d+)\s+data\s+follows\s*$/)
      {     $nbline = $1;
            if($phase!=2)
            {   print "error while reading data\n";
                exit;
            }
            for (my $i=0;$i<$nbline;$i++)
            {     my $line=<DF>;
                  chomp($line);
                  my @temp=split(/\s+/,$line);
                  shift @temp if (scalar(@temp)>8);
                  $hex[$numhex] = \@temp;
                  $numhex++;
            }
            $phase=3;
      }
      else
      {     if($phase!=3)
            {   print "error while reading data\n";
                exit;
            }
      }
}
close (DF);
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

my ($xmin,$xmax,$ymin,$ymax,$zmin,$zmax) = ($hugeval,-$hugeval,$hugeval,-$hugeval,$hugeval,-$hugeval);
foreach (@nodes)
{   $xmin=$_->[0] if ($_->[0] < $xmin);
    $ymin=$_->[1] if ($_->[1] < $ymin);
    $zmin=$_->[2] if ($_->[2] < $zmin);
    $xmax=$_->[0] if ($_->[0] > $xmax);
    $ymax=$_->[1] if ($_->[1] > $ymax);
    $zmax=$_->[2] if ($_->[2] > $zmax);
}

print "the formula to compute the total number of GLL points in the two-level superbrick is the following :\n";
print "number of elements : $numhex\n";
my $faces = get_connexion_faces(\@cubes);
$nbfaces = scalar(@$faces);
print "number of contact surfaces between these elements : ",$nbfaces,"\n";
my $edges = get_connexion_edges($faces);
#$faces=undef;
$edges = edge_clean_up($edges);
$nbedges = scalar(@$edges);
print "number of (internal) contact edges between these surfaces : ",$nbedges,"\n";
my $points = get_connexion_points($edges);
#$edges=undef;
$points = point_clean_up($points);
$nbpoints = scalar(@$points);
#$points=undef;
print "number of (internal) contact points between these edges : ",$nbpoints,"\n";
print "\ntherefore the final formula is\nNGLOB = $numhex.NGLL^3 - $nbfaces.NGLL^2 + $nbedges.NGLL - $nbpoints\n";
$res = (($numhex*($NGLL**3))-($nbfaces*($NGLL**2))+($nbedges*$NGLL)-$nbpoints);
print ("for NGLL = $NGLL, NGLOB = $res\n") if ($NGLL);
exit;

sub point_clean_up      # remove the points in double
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

sub edge_clean_up       # remove the edges in double
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

sub get_connexion_points # in an edges set
{     my $edges = shift;
      my @conn_points;
      my $nb=0;
      foreach $couple (@$edges)
      {     $nb++;
            for (my $i=$nb;$i<$nbedges;$i++)
            {     my $ok = 0;
                  my $tmp = get_conn_point($couple, $edges->[$i]);
                  push (@conn_points, $tmp) if ($tmp);
            }
      }
      return \@conn_points;
}

sub get_conn_point  # between 2 edges
{     my ($couple_1, $couple_2) = @_;
      foreach $pt1 (@$couple_1)
      {     foreach $pt2 (@$couple_2)
            {     return $pt1 if (($pt1->[0] == $pt2->[0] and $pt1->[1] == $pt2->[1] and $pt1->[2] == $pt2->[2]) and (is_in($pt1)));
            }
      }
      return 0;
}

sub get_connexion_edges  # in an faces set
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

sub get_face_edges  # get the 4 edges of a quad
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

sub get_conn_edges  # between 2 faces
{     my ($edges_set_1, $edges_set_2) = @_;
      foreach $edge1 (@$edges_set_1)
      {     foreach $edge2 (@$edges_set_2)
            {     if (eqedge($edge1,$edge2))
                  {     return $edge1 unless (is_external($edge1));
                  }
            }
      }
      return 0;
}

sub is_external   # check if an edge is at the surface of the brick
{   my $edge = shift;
    my($x1,$x2,$y1,$y2,$z1,$z2) = ($edge->[0][0],$edge->[1][0],$edge->[0][1],$edge->[1][1],$edge->[0][2],$edge->[1][2]);
    return 1 if(($x1==$x2 and ($x1==$xmin or $x1==$xmax)) or ($y1==$y2 and ($y1==$ymin or $y1==$ymax)) or ($z1==$z2 and ($z1==$zmin or $z1==$zmax)));
    return 0;
}

sub is_in     # check if a point is inside the brick, not at its surface
{   my $point = shift;
    my ($x,$y,$z) = ($point->[0],$point->[1],$point->[2]);
    return 0 if ($x==$xmin or $x==$xmax or $y==$ymin or $y==$ymax or $z==$zmin or $z==$zmax);
    return 1;
}

sub eqedge    # check if 2 edges are equal or not
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

sub get_connexion_faces # for an hexahedron set
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

sub get_cube_faces # get the 6 faces of an hexahedron
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

sub get_conn_faces # between 2 hexahedron
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

sub eqface # check if 2 faces are equal or not
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
