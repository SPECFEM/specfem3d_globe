#! /usr/bin/perl

# give the modules & subroutines dependencies for each file of the SPECFEM3D_GLOBE code.
#
# regular expressions may have to be tuned.
#
# David Michea, University of Pau, France, december 2006


$|++;

if (scalar(@ARGV)<2)
{ print "usage : $0 fileout source_dir\n";
  exit;
}
my $fileout = $ARGV[0];
my $src_dir = $ARGV[1];
my @mystruct;
print ">> data structure creation\n";
my @files = glob("$src_dir/*.f90");
my $nbfiles = scalar(@files);
my $cpt=0;
foreach my $fichier (@files)
{ $cpt++;
  my $line_num = 0;
  my ($pfiletab, $pmods, $psubs, $pmodst, $psubst) = ([],[],[],[],[]);
  my %hmod = ();
  my %hsub = ();
  my %hmodt = ();
  my %hsubt = ();
  open (DF, $fichier);
  while (<DF>)
  { $line_num++;
    if ($_ =~ m/^\s*use\s+([a-zA-Z_0-9-]*)$/) {
      push (@$pmods, [$1]) unless(exists($hmod{$1}));
      $hmod{$1}++;
    }
    elsif (($_ =~ m/^(?:[^!]*\s+|\s*)call\s+([a-zA-Z_0-9-]*).*$/) && ($_ !~ m/MPI_[A-Z]+/)) {
      push (@$psubs, [$1]) unless(exists($hsub{$1}));
      $hsub{$1}++;
    }
    elsif ($_ =~ m/^\s*module\s+([a-zA-Z_0-9-]*).*$/) {
      push (@$pmodst, [$1, $line_num]) unless(exists($hmodt{$1}));
      $hmodt{$1}++;
    }
    elsif ($_ =~ m/^\s*subroutine\s+([a-zA-Z_0-9-]*).*$/) {
      push (@$psubst, [$1, $line_num]) unless(exists($hsubt{$1}));
      $hsubt{$1}++;
    }
  }
  close (DF);
  $fichier =~ s#.*/(.*)$#$1#g;
  push (@$pfiletab, $fichier, $pmods, $psubs, $pmodst, $psubst);
  push (@mystruct, $pfiletab);
  print "\r".int(($cpt/$nbfiles)*100)." %";
}
print "\n>> analysis\n";
$cpt=0;
foreach my $ptrfile (@mystruct)
{ $cpt++;
  print "\r".int(($cpt/$nbfiles)*100)." %";
  for (my $i=0;$i<2;$i++)
  { foreach my $ptrmod (@{$ptrfile->[1+$i]})
    { my $last = 0;
      foreach my $ptrfilet (@mystruct)
      { foreach my $ptrmodt (@{$ptrfilet->[3+$i]})
        { if ($ptrmod->[0] eq $ptrmodt->[0])
          { push (@$ptrmod, $ptrfilet->[0]." [l ".$ptrmodt->[1]."]");
            $last++;
            last;
          }
        }
        last if ($last>0);
      }
    }
  }

}
print "\n>> results writing in $fileout\n";
open (OUT, ">$fileout");
$cpt=0;
foreach my $ptrfile (@mystruct)
{ $cpt++;
  if (scalar(@{$ptrfile->[1]}) > 0 || scalar(@{$ptrfile->[2]}) > 0) {
    my $l = int((76 - length($ptrfile->[0]))/2);
    print OUT "x"x$l,"- ".$ptrfile->[0]." -","x"x($l+((76 - length($ptrfile->[0]))%2)),"\n";
    print OUT ">> MODULES\n" if (scalar(@{$ptrfile->[1]}));
    foreach my $ptrmod (@{$ptrfile->[1]})
    { print OUT "\t".$ptrmod->[0]."  ->  ".$ptrmod->[1]."\n";
    }
    print OUT ">> SUBROUTINES\n" if (scalar(@{$ptrfile->[2]}));
    foreach my $ptrsub (@{$ptrfile->[2]})
    { print OUT "\t".$ptrsub->[0]."  ->  ".$ptrsub->[1]."\n";
    }
    print OUT "\n\n";
  }
  print "\r".int(($cpt/$nbfiles)*100)." %";
}
print "\n\ndone\n";
close(OUT);
exit;
