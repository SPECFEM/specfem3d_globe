#!/usr/bin/perl

use Getopt::Std;
use DBI;
use strict;

$|++;
# for using database analysis, you have to install, run & configure mysql-server

# constantes utilisees globalement
my $pi = 3.141592653589793;
my $pisur2 = (3.141592653589793 / 2);
my $pisur4 = (3.141592653589793 / 4);
my $Ricb = 1221;
my $nbnode=8;
my $hugeval = 1000000;
my $NEXXI = 256;
my $nb_class_histo = 20;
my @iaddx = (0,1,2,2,2,1,0,0);
my @iaddy = (0,0,0,1,2,2,2,1);
my @order=([0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,0]);
my @mestitres =('taille + grde arete\n-------------------\ntaille + petite arete',"aspect ratio : pire element","aspect ratio moyen","heuristique pire skewness","heuristique skewness moy","moyenne des pires skewness","moyenne des skewness moy","pire des pires skewness","pire des skewness moy");
my $dbh;
my $database = "dbi:mysql:MY_DB:127.0.0.1";
my $dbUser = "root";
my $dbPasse = "mynewpassword";
# ---------------

# <main>
my %opts;
getopts("avdc", \%opts) or print_usage();
if (exists($opts{c}))
{ print_usage() if (scalar(@ARGV)!=1);
  my $max_angle_param = shift(@ARGV);
  if ($max_angle_param<90 or $max_angle_param>180)
  { print "angle max have to belong to [90 degrees,180 degrees]\n";
    print_usage();
  }
  my $alpha=1;
  while ((my $angle_max = (($pisur2 + 2 * atan2($alpha,1))/(2*$pi))*360)> $max_angle_param)
  { $alpha-=0.00001;}
  print "alpha=$alpha\n";
  exit;
}
if (exists($opts{v}) and not exists($opts{a}))
{ print_usage() if (scalar(@ARGV)!=3);
  my ($Rcube,$alpha,$filename) = @ARGV;
  chomp($filename);
  if ($alpha<0 or $alpha>1)
  { print "alpha have to belong to [0,1]\n";
    print_usage();
  }
  if ($Rcube>$Ricb)
  { print "Rcube have to be < Ricb=$Ricb\n";
    print_usage();
  }
  visualize_mesh($Rcube,$alpha,$filename);
  print "done.\n";
}
elsif (exists($opts{a}) and not exists($opts{v}))
{ print_usage() if (scalar(@ARGV)!=7);
  my ($filename,$alpha_min,$alpha_max,$alpha_step,$Rcube_min,$Rcube_max,$Rcube_step) = @ARGV;
  chomp($Rcube_step);
  $dbh = DBI->connect($database, $dbUser, $dbPasse) if (exists($opts{d}));
  #db_create();
  analyse($filename,$alpha_min,$alpha_max,$alpha_step,$Rcube_min,$Rcube_max,$Rcube_step);
  $dbh->disconnect if (exists($opts{d}));
  print "done.\n";
}
else
{ print_usage();
}
# </main>

sub print_usage
{ print "\nusage\n\noption -v : visualize a mesh for a radius of central cube (km) and alpha [0,1]\n";
  print "option -a : analyse skewness & aspect ratio for all meshes\n";
  print "option -c : give alpha for a given max angle\n\n";
  print "$0\t-v R_central_cube alpha fileout.dx\n";
  print "\t|\t-c max_angle\n";
  print "\t|\t-a fileout alpha_min alpha_max alpha_step Rcube_min Rcube_max Rcube_step\n\n";
  exit;
}

sub visualize_mesh
{ my ($Rcube,$alpha,$filename) = @_;
  my $nex_xi = $NEXXI/16;
  my $nex_eta = compute_ner($nex_xi,$Ricb,$Rcube,$alpha);
  my $nspec_cube=0;
  my $nspec_chunks=0;
  my $points = mesh($Rcube,$alpha,$nex_xi,$nex_eta,\$nspec_cube, \$nspec_chunks);
  writedxfile($points,"$filename",$nex_xi,$nex_eta,$nspec_cube,$nspec_chunks);
}

sub analyse
{ my ($fileout,$alpha_min,$alpha_max,$alpha_step,$Rcube_min,$Rcube_max,$Rcube_step) = @_;
  open(DATA, ">$fileout".".dat");
  print "progression :\n";
  #exit;
  for (my $ialpha=$alpha_min;$ialpha<$alpha_max;$ialpha+=$alpha_step)
  { for (my $iRcube=$Rcube_min;$iRcube<$Rcube_max;$iRcube+=$Rcube_step)
    { # variables  pour l'analyse
      my $max_edge_of_mesh=-1;
      my $min_edge_of_mesh = $hugeval;
      my $aspect_ratio_max = -1;
      my $sum_aspect_ratio = 0;
      my @skewness_worth_histo;
      my @skewness_moy_histo;
      my $skewness_w_w;
      my $skewness_w_m;
      my $skewness_m_w;
      my $skewness_m_m;

      my $nex_xi = $NEXXI/16;
      # calcul du nb d'elements radiaux dans les chunks
      my $nex_eta = compute_ner($nex_xi,$Ricb,$iRcube,$ialpha);
      my $nspec_cube=0;
      my $nspec_chunks=0;
      # maillage
      my $points = mesh($iRcube,$ialpha,$nex_xi,$nex_eta,\$nspec_cube, \$nspec_chunks);
      # analyse des elements
      for (my $ispec=1;$ispec<=($nspec_cube+$nspec_chunks);$ispec++)
      { my @elem=get_elem($points,$ispec);
        my ($edgemin,$edgemax) = get_size_min_max(@elem);
        my $aspect_ratio = $edgemax/$edgemin;
        my ($skewness_moy,$skewness_worth) = get_skewness_moy_worth(get_angles(@elem));

        $max_edge_of_mesh=$edgemax if ($edgemax>$max_edge_of_mesh);
        $min_edge_of_mesh=$edgemin if ($edgemin<$min_edge_of_mesh);
        $aspect_ratio_max = $aspect_ratio if ($aspect_ratio>$aspect_ratio_max);
        $sum_aspect_ratio+=$aspect_ratio;
        $skewness_worth_histo[get_histo_class($skewness_worth)]++;
        $skewness_moy_histo[get_histo_class($skewness_moy)]++;

        $skewness_w_w = $skewness_worth if ($skewness_worth>$skewness_w_w);
        $skewness_m_w += $skewness_worth;

        $skewness_w_m= $skewness_moy if ($skewness_moy>$skewness_w_m);;
        $skewness_m_m += $skewness_moy;
      }
      $skewness_m_w/=($nspec_cube+$nspec_chunks);
      $skewness_m_m/=($nspec_cube+$nspec_chunks);
      # les trois aspects de l'analyse
      my $f1=$max_edge_of_mesh/$min_edge_of_mesh;
      my @f2=($aspect_ratio_max,($sum_aspect_ratio/($nspec_cube+$nspec_chunks)));
      my @f3=(\@skewness_worth_histo,\@skewness_moy_histo);

      # ecriture des resultats
      if (exists($opts{d}))
      { #db_insert();
      }
      else
      { print DATA join(" ",($ialpha,$iRcube,$f1,join(" ",@f2),repart_heurist($f3[0]),repart_heurist($f3[1])),$skewness_m_w,$skewness_m_m,$skewness_w_w,$skewness_w_m),"\n";
      }
      print ("\r", int((($iRcube-$Rcube_min)/$Rcube_step)+((($ialpha-$alpha_min)/$alpha_step)*(($Rcube_max-$Rcube_min)/$Rcube_step))+1)," / ",((($alpha_max-$alpha_min)/$alpha_step)*(($Rcube_max-$Rcube_min)/$Rcube_step)));
    }
    print DATA "\n";
  }
  print "\ndone.\nlaunching gnuplot\n";
  write_gnuplot_prog($fileout,($alpha_max-$alpha_min)/$alpha_step,($Rcube_max-$Rcube_min)/$Rcube_step,@mestitres);
  close(DATA);
  #create_dx_visu($fileout,$alpha_min,$alpha_max,$alpha_step,$Rcube_min,$Rcube_max,$Rcube_step,0,3);
  exec("gnuplot $fileout.gplot");
# create_dx_visu($fileout,$min_alpha,$max_alpha,$min_Rcube,$max_Rcube,$step_alpha,$step_Rcube,$nbclass);
}

sub write_gnuplot_prog
{ my ($filebase,$ech_alpha,$ech_Rcube,@tittles) = @_;
  open(OUT,">$filebase".".gplot");
  my $nbcol=scalar(@tittles);
  for (my $icol=3;$icol<($nbcol+3);$icol++)
  { print OUT "set hidden3d\n";
    print OUT "set isosamples $ech_alpha,$ech_Rcube\n";
    print OUT "set grid\n";
    print OUT "set title \"",$tittles[$icol-3],"\"\n";
    print OUT "set data style line\n";
    print OUT "set contour #both\n";
    print OUT "splot '$filebase.dat' using 1:2:$icol\n";
    print OUT "pause -1\n";
  }
  for (my $icol=3;$icol<($nbcol+3);$icol++)
  { print OUT "set hidden3d\n";
    print OUT "set isosamples $ech_alpha,$ech_Rcube\n";
    print OUT "set grid\n";
    print OUT "set title \"",$tittles[$icol-3],"\"\n";
    print OUT "set contour\n";
    print OUT "set cntrparam levels 20\n";
    print OUT "set data style line\n";
    print OUT "set nosurface\n";
    print OUT "set view 0,0\n";
    print OUT "splot '$filebase.dat' using 1:2:$icol\n";
    print OUT "pause -1\n";
  }
  close(OUT);
  return;
}


sub repart_heurist
{ my $histo=shift;
  my $ret;
  my $cpt=$nb_class_histo;
  my $nbelem=0;
  foreach (@$histo)
  { $nbelem+=$_;
  }
  foreach (@$histo)
  { $ret+=($_/$nbelem)*$cpt;
    $cpt--;
  }
  return $ret;
}

sub get_histo_class
{ my $skew=shift;
  for (my $class=0;$class<$nb_class_histo;$class++)
  { return $class if ($skew<(($class+1)*(1/$nb_class_histo)));
  }
}

sub get_skewness_moy_worth
{ my (@angles) = @_;
  my $somme=0;
  my $max_skew=-1;
  foreach (@angles)
  { my $cur_skew = abs((2*$_-$pi)/$pi);
    $somme+=$cur_skew;
    $max_skew=$cur_skew if ($cur_skew>$max_skew);
  }
  return (($somme/scalar(@angles)),$max_skew);
}

sub get_elem
{ my ($points,$ispec) = @_;
  my @element = @$points[(($ispec-1)*$nbnode)..($ispec*$nbnode-1)];
  return @element;
}

sub get_angles
{   my (@elem) = @_;
  my @angles;
  push(@elem,(@elem[0,1]));
  for (my $iedge=0;$iedge<$nbnode/2;$iedge++)
  { my ($point2,$point3,$point4) = @elem[($iedge*2+1)..($iedge*2+3)];
    my($x2,$y2,$x3,$y3,$x4,$y4) = (split(/ /,$point2),split(/ /,$point3),split(/ /,$point4));
    my ($xv1,$yv1,$xv2,$yv2) = ($x2-$x3,$y2-$y3,$x4-$x3,$y4-$y3);
    my $angle = acos(($xv1*$xv2+$yv1*$yv2)/(sqrt($xv1**2+$yv1**2) * sqrt($xv2**2+$yv2**2)));
    push(@angles,$angle);
  }
  return @angles;
}

sub acos {
   my($x) = @_;
   my $ret = atan2(sqrt(1 - $x**2), $x);
   return $ret;
}

sub get_size_min_max
{ my (@elem) = @_;
  push(@elem,$elem[0]);
  my $sizemax=-1;
  my $sizemin=$hugeval;
  for (my $iedge=0;$iedge<$nbnode/2;$iedge++)
  { my ($point1,$point2,$point3) = @elem[($iedge*2)..($iedge*2+2)];
    my($x1,$y1,$x2,$y2,$x3,$y3) = (split(/ /,$point1),split(/ /,$point2),split(/ /,$point3));
    my $size= dist($x1,$y1,$x2,$y2) + dist($x2,$y2,$x3,$y3);
    $sizemax=$size if ($size>$sizemax);
    $sizemin=$size if ($size<$sizemin);
  }
  return ($sizemin,$sizemax);
}

sub dist
{ my ($x1,$y1,$x2,$y2) = @_;
  my $res=sqrt(($x1-$x2)**2+($y1-$y2)**2);
  if ($res==0)
  { print "\n$x1 , $y1 , $x2 , $y2\n";
    exit;
  }
  return $res;

}

sub mesh
{ my ($Rcube,$alpha,$nex_xi,$nex_eta,$nspec_cube, $nspec_chunks) = @_;
  my @points;
  # maillage des 4 chunks
  for (my $ichunk=0;$ichunk<4;$ichunk++)
  { for (my $ix=0;$ix<$nex_xi*2;$ix+=2) # boucle circulaire
    { for (my $iy=0;$iy<$nex_eta*2;$iy+=2) # boucle radiale
      { for (my $inode=0;$inode<$nbnode;$inode++)
        { my ($x,$y) = compute_coord_chunk($ix+$iaddx[$inode], $iy+$iaddy[$inode],$nex_xi*2,$nex_eta*2,$Rcube,$ichunk,$alpha);
          push (@points, join(" ",($x,$y)));
        }
        $$nspec_chunks++;
      }
    }
  }
  # maillage du carre central
  for (my $ix=0;$ix<2*$nex_xi;$ix+=2)
  { for (my $iy=0;$iy<2*$nex_xi;$iy+=2)
    { for (my $inode=0;$inode<$nbnode;$inode++)
      { my ($x,$y) = compute_coord_central_cube($ix+$iaddx[$inode], $iy+$iaddy[$inode],$nex_xi*2,$nex_xi*2,$Rcube,$alpha);
        push (@points, join(" ",($x,$y)));
      }
      $$nspec_cube++;
    }
  }
  return \@points;
}

sub compute_coord_chunk
{ my ($ix,$iy,$nb_elem_x,$nb_elem_y,$Rcube,$ichunk,$alpha) = @_;
  my $ratio_x = $ix/$nb_elem_x;
  my $ratio_y = $iy/$nb_elem_y;

  my $factx = (2*$ratio_x -1);
  my $xi = $pisur2*$factx;

  # coordonnees des extremites d'aretes a la surface du CC
  my $xcc=($Rcube / sqrt(2)) * $factx;
  my $ycc=($Rcube / sqrt(2)) * (1 + cos($xi)*$alpha / $pisur2);
  # coordonnees des extremites d'aretes a la surface de l'ICB
  my $xsurf = $Ricb * cos(3*$pisur4 - $ratio_x * $pisur2);
  my $ysurf = $Ricb * sin(3*$pisur4 - $ratio_x * $pisur2);

  my $deltax=$xsurf-$xcc;
  my $deltay=$ysurf-$ycc;
  # coordonnees du point
  my $x = $xsurf-$ratio_y*$deltax;
  my $y = $ysurf-$ratio_y*$deltay;
  # changement de repere pour les trois autres chunks
  if ($ichunk==1)
  { my $temp=$x; $x=$y; $y=-$temp;
  }
  elsif ($ichunk==2)
  { $x=-$x; $y=-$y;
  }
  elsif ($ichunk==3)
  { my $temp=$x; $x=-$y; $y=$temp;
  }
  return ($x,$y);
}

sub compute_coord_central_cube
{ my ($ix,$iy,$nb_elem_x,$nb_elem_y,$radius,$alpha) = @_;
  my ($x,$y);

  my $ratio_x = $ix/$nb_elem_x;
  my $ratio_y = $iy/$nb_elem_y;

  my $factx = (2*$ratio_x -1);
  my $facty = (2*$ratio_y -1);

  my $xi = $pisur2*$factx;
  my $eta = $pisur2*$facty;

  $x=($radius / sqrt(2)) * $factx  * (1 + cos($eta)*$alpha / $pisur2);
        $y=($radius / sqrt(2)) * $facty  * (1 + cos($xi)*$alpha / $pisur2);

  return ($x,$y);
}

sub compute_ner
{ my ($nex_xi,$Ricb,$Rcube,$alpha) = @_;

  my $somme=0;
  # on fait la moyenne des distances ICB<>CC sur 1/2 chunk (suffisant)
  # attention ! valable seulement si $nex_xi est pair, il faudra prendre en compte le cas general si besoin
  for (my $ix=0;$ix<=$nex_xi/2;$ix++)
  {
    my $ratio_x = $ix/$nex_xi;
    my $factx = 2*$ratio_x-1;
    my $xi = $pisur2*$factx;

    # coordonnees des extremites d'aretes a la surface du CC
    my $x=($Rcube / sqrt(2)) * $factx;
    my $y=($Rcube / sqrt(2)) * (1 + cos($xi)*$alpha / $pisur2);
    # coordonnees des extremites d'aretes a la surface de l'ICB
    my $xsurf = $Ricb * cos(3*$pisur4 - $ratio_x * $pisur2);
    my $ysurf = $Ricb * sin(3*$pisur4 - $ratio_x * $pisur2);
    # taille de l'arete
    my $dist_cc_icb = sqrt(($xsurf-$x)**2+($ysurf-$y)**2);
    # on double le poids de chaque arete sauf l'arete centrale (la derniere) a ne compter qu'une fois
    $dist_cc_icb*=2 unless($ix==$nex_xi/2);
    $somme+=$dist_cc_icb;
  }
  my $dist_moy = $somme/($nex_xi+1);
  my $ner = arrondi($dist_moy/(($pi*$Ricb)/(2*$nex_xi)));
  return $ner;
}

sub arrondi
{ my $input = shift;
  my $arrondi = sprintf("%.0f", $input);
  return $arrondi;
}

sub writedxfile
{ my ($points,$filename,$nb_elem_x,$nb_elem_y,$nspec_cube,$nspec_chunks) = @_;

  my $nblignes = scalar(@$points);
  open(OUT, ">$filename");
  print OUT " object 1 class array type float rank 1 shape 2 items         ",$nblignes,"  data follows\n";
  foreach (@$points)
  { print OUT $_,"\n";
  }
  print OUT " object 2 class array type int rank 1 shape 2 items         ",$nblignes,"  data follows\n";
  for (my $i=0;$i<($nspec_cube+$nspec_chunks);$i++)
  { for (my $j=0;$j<$nbnode;$j++)
    { print OUT ($i*$nbnode+$order[$j]->[0]," ",$i*$nbnode+$order[$j]->[1],"\n");
    }
  }
  print OUT <<"EOF"
  attribute \"element type\" string \"lines\"
 attribute \"ref\" string \"positions\"
 object 3 class array type float rank 0 items $nblignes  data follows
EOF
;
  my $cpt=0;
  foreach (@$points)
  { $cpt++;
    if ($cpt<=($nspec_chunks*$nbnode))
    { print OUT "1\n" ;}
    else
    { print OUT "2\n" ;}
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
}

sub create_dx_visu
{ my ($filein,$min_alpha,$max_alpha,$step_alpha,$min_Rcube,$max_Rcube,$step_Rcube,$nbclass,$indice) = @_;

  open(IN,"$filein".".dat");
  my $nblignes = (($max_alpha-$min_alpha)/$step_alpha)*(($max_Rcube-$min_Rcube)/$step_Rcube);

  open(OUT, ">$filein".".dx");
  print OUT " object 1 class array type float rank 1 shape 2 items         ",$nblignes,"  data follows\n";
  for (my $ialpha=$min_alpha;$ialpha<$max_alpha;$ialpha+=$step_alpha)
  { for (my $iRcube=$min_Rcube;$iRcube<$max_Rcube;$iRcube+=$step_Rcube)
    { print OUT "$ialpha $iRcube\n";
    }
  }
  print OUT " object 2 class array type int rank 1 shape 2 items         ",$nblignes,"  data follows\n";
  for (my $ialpha=0;$ialpha<($max_alpha-$min_alpha);$ialpha+=$step_alpha)
  { for (my $iRcube=0;$iRcube<($max_Rcube-$min_Rcube);$iRcube+=$step_Rcube)
    { print OUT $ialpha+($max_alpha-$min_alpha)*$iRcube," ";
      print OUT ($ialpha+1)+($max_alpha-$min_alpha)*$iRcube," ";
      print OUT $ialpha+($max_alpha-$min_alpha)*($iRcube+1)," ";
      print OUT ($ialpha+1)+($max_alpha-$min_alpha)*($iRcube+1),"\n";
    }
  }
  print OUT <<"EOF"
  attribute \"element type\" string \"squares\"
 attribute \"ref\" string \"positions\"
 object 3 class array type float rank 0 items $nblignes  data follows
EOF
;
  while(my $ligne = <IN>)
  { my @tokens = split(";",$ligne);
    my $zvalue = $tokens[$indice];
    print OUT int($zvalue*100),"\n";
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
}

