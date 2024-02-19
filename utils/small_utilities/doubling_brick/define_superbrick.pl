#! /usr/bin/perl

use Data::Dumper;

#<main>
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
my @nodes;
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

my $compt=0;
foreach (@nodes)
{   $compt++;
    print "x_superbrick($compt) = ",$_->[0],"\n";
    print "y_superbrick($compt) = ",$_->[1],"\n";
    print "z_superbrick($compt) = ",$_->[2],"\n";
    print "\n";
}
reord_points(\@cubes);
my @bound = get_boundaries(\@cubes);
my $compt_elem==0;
print "iboun_sb(:,:) = .false.\n";
foreach my $res (@bound)
{ $compt_elem++;
  my $compt_boun==0;
  foreach (@$res)
  { $compt_boun++;
    print "iboun_sb($compt_elem,$compt_boun) = .true.\n" if ($_ == 1);
  }
}
exit;
#</main>


sub get_boundaries
{ my $cubes = shift;
  my $nb=0;
  my @bound;
  foreach $hex (@$cubes)
  { $nb++;
    my $cur_cube_faces = get_cube_faces($hex);
    my @res =get_pos_faces($cur_cube_faces);
    push(@bound,\@res);
  }
  return @bound;
}

sub get_pos_faces
{ my $faces = shift;
  my ($xmin,$xmax,$ymin,$ymax,$bottom,$top) = (0,0,0,0,0,0);
  foreach my $face (@$faces)
  { my ($moy_x,$moy_y,$moy_z) = analyse_face2($face);
    $xmin=1 if ($moy_x == 0);
    $xmax=1 if ($moy_x == 2);
    $ymin=1 if ($moy_y == 0);
    $ymax=1 if ($moy_y == 2);
    $bottom=1 if ($moy_z == 0);
    $top=1 if ($moy_z == 2);
  }
  return ($xmin,$xmax,$ymin,$ymax,$bottom,$top);
}

sub analyse_face2
{   my $face = shift;
    my $moy_x=0;
    my $moy_y=0;
    my $moy_z=0;
    foreach my $corner (@$face)
    {   $moy_x+=$corner->[0];
        $moy_y+=$corner->[1];
        $moy_z+=$corner->[2];
    }
    $moy_x/=4;
    $moy_y/=4;
    $moy_z/=4;
    return ($moy_x,$moy_y,$moy_z);
}


sub reord_points
{     my $cubes = shift;
      my @conn_faces;
      my $nb=0;
      foreach $hex (@$cubes)
      {     $nb++;
            my $cur_cube_faces = get_cube_faces($hex);
            my $facebas = get_bas($cur_cube_faces);
            my $facehaut = get_haut($cur_cube_faces);
            my @points_ord = get_ord ($facebas, $facehaut);
            for (my $i=1;$i<9;$i++)
            {   print "ibool_superbrick($i, $nb) = ",shift(@points_ord),"\n";
            }
            print "\n";
      }
      return
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

sub get_ord
{   my ($faceb, $faceh) = @_;
    my @bas = ordonne($faceb);
    my @haut = ordonne($faceh);
    my @ord = (@bas,@haut);
    return @ord;
}

sub ordonne
{   my $face = shift;
    my @ord;
    my $xmin = 1000;
    my $xmax = -1000;
    my $ymin = 1000;
    my $ymax = -1000;
    foreach my $corner (@$face)
    {   $xmin = $corner->[0] if ($corner->[0] < $xmin);
        $ymin = $corner->[1] if ($corner->[1] < $ymin);
        $xmax = $corner->[0] if ($corner->[0] > $xmax);
        $ymax = $corner->[1] if ($corner->[1] > $ymax);
    }
    foreach my $corner (@$face)
    {   $ord[0] = get_num_point($corner) if ($corner->[0] == $xmin and $corner->[1] == $ymin);
        $ord[1] = get_num_point($corner) if ($corner->[0] == $xmax and $corner->[1] == $ymin);
        $ord[2] = get_num_point($corner) if ($corner->[0] == $xmax and $corner->[1] == $ymax);
        $ord[3] = get_num_point($corner) if ($corner->[0] == $xmin and $corner->[1] == $ymax);
    }
    return @ord;
}

sub get_num_point
{   my $point = shift;
    my $x = $point->[0];
    my $y = $point->[1];
    my $z = $point->[2];
    for (my $i=0;$i<$numnode;$i++)
    {   return ($i+1) if ($x == $nodes[$i][0] and $y == $nodes[$i][1] and $z == $nodes[$i][2]);
    }
    print "point inconnu ! $x, $y, $z\n";
    exit;
}

sub analyse_face
{   my $face = shift;
    my $vert=0;
    my $moy_x=0;
    my $moy_y=0;
    my $moy_z=0;
    foreach my $corner (@$face)
    {   $moy_x+=$corner->[0];
        $moy_y+=$corner->[1];
        $moy_z+=$corner->[2];
    }
    $moy_x/=4;
    $moy_y/=4;
    $moy_z/=4;
    $vert=1 if (($moy_x == $face->[0][0]) or ($moy_y == $face->[0][1]));
    return ($moy_z,$vert)
}

sub get_bas
{   my $faces = shift;
    my $min_z = 1000;
    my $numfacemin = 100;
    my $cpt = 0;
    foreach my $face (@$faces)
    {   my ($moy_z, $vert) = analyse_face($face);
        if (($moy_z < $min_z) or (($moy_z == $min_z) and (not $vert)))
        {   $min_z = $moy_z;
            $numfacemin = $cpt;
        }
        $cpt++;
    }
    return $faces->[$numfacemin];
}

sub get_haut
{   my $faces = shift;
    my $max_z = -1000;
    my $numfacemax = 100;
    my $cpt = 0;
    foreach my $face (@$faces)
    {   my ($moy_z, $vert) = analyse_face($face);
        if (($moy_z > $max_z) or (($moy_z == $max_z) and (not $vert)))
        {   $max_z = $moy_z;
            $numfacemax = $cpt;
        }
        $cpt++;
    }
    return $faces->[$numfacemax];
}
