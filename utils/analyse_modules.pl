#! /usr/bin/perl

# search the incoherences in parameters passing or type declaration for the new structures created in place of the modules in the SPECFEM3D_GLOBE code.
# structures are declared with keyword "type"
# no output if no problem.

# this script has been written for this particular usage, not for a general purpose.

# David Michea. University of Pau. February 2007.

$|++;

if (scalar(@ARGV)<1)
{ print "usage : $0 source_dir\n";
  exit;
}
my $src_dir = $ARGV[0];
chomp($src_dir);
if ($src_dir =~m /^(.*)\/$/)
{ $src_dir =~ $1;
}
my @mystruct;
my %htypes;
print ">> data structure creation\n";
my @files = glob("$src_dir/*.f90");
my $nbfiles = scalar(@files);
my $cpt=0;
foreach my $fichier (@files)  # parcours de chaque fichier
{ $cpt++;
  my $line_num = 0;
  my $pfiletab = [];
  open (DF, $fichier);
  $fichier =~ s#.*/(.*)$#$1#g;  #
  push(@$pfiletab,$fichier);
  while (<DF>)    # parcours de chaque ligne
  { $line_num++;
  # analyse de chaque declaration de subroutine
    # declaration sur une ligne
    if ($_ =~ m/^\s*subroutine\s+([a-zA-Z_0-9-]*)\s*(?:\(([^!]*)\))?(?:!.*|\s*)$/) {
      my $nom = $1;
      my @params = split(/,/,$2);
      @params = map(trim($_),@params);
      push (@$pfiletab,analyse_sub($fichier,$nom,\$line_num,@params));
    }
    # declaration sur plusieurs lignes (&)
    elsif ($_ =~ m/^\s*subroutine\s+([a-zA-Z_0-9-]*)\s*\((.*)&\s*$/)
    { my $nom = $1;
                        my @params = split(/,/,$2);
      pop(@params) if ($params[-1] =~ m/^\s*$/);
      while(<DF>)
      { $line_num++;
        if ($_ =~ m/^\s*(.*)&\s*$/)
        { push(@params,split(/,/,$1));
          pop(@params) if ($params[-1] =~ m/^\s*$/);
        }
        elsif ($_ =~ m/^\s*(.*)\)\s*$/)
        { push(@params,split(/,/,$1));
                                        pop(@params) if ($params[-1] =~ m/^\s*$/);
          last;
        }
        else
        { print "probleme : declaration non conforme dans $fichier l $line_num\n";exit;
        }
      }
                        @params = map(trim($_),@params);
      push (@$pfiletab,analyse_sub($fichier,$nom,\$line_num,@params));
    }
    else
    {#  print "$_";
    }
  }
  close (DF);
  push (@mystruct, $pfiletab);
  print "\r".int(($cpt/$nbfiles)*100)." %";
}
print("\n");
#exit;

print "\n>> structure analysis\n";
foreach my $ptrfile (@mystruct)
{       my $file = @$ptrfile[0];
  my $i=0;
  foreach my $ptrsub (@$ptrfile)
  {
    $i++;
    next if ($i==1);
    my $nomsub = $ptrsub->[0];
    my $ptrparamsin = $ptrsub->[1];
    my $ptrtypesdec = $ptrsub->[2];
    test_paramsin($ptrparamsin,$ptrtypesdec,$nomsub,$file);
    foreach my $ptrfile2  (@mystruct)
    { my $file2 = @$ptrfile2[0];
      my $j=0;
      foreach my $ptrsub2 (@$ptrfile2)
      {
        $j++;
        next if ($j==1);
        my $nomsub2 = $ptrsub2->[0];
                    my $ptrparamsin2 = $ptrsub2->[1];
        if ($#{$ptrsub2}>3)
        { foreach $ptrsubsub (@$ptrsub2[3,$#{$ptrsub2}])
          {
            my $nomsubsub = $ptrsubsub->[0];
            my $lignesubsub = $ptrsubsub->[1];
                        my $ptrparamsinsubsub = $ptrsubsub->[2];
            if ($nomsubsub eq $nomsub)
            { if (scalar(@$ptrparamsin) != scalar(@$ptrparamsinsubsub))
              { print "for subroutine $nomsub, declared in file $file, its call in the subroutine $nomsub2 line $lignesubsub of file $file2";
                print " is incomplete :\n$nomsub(",join(",",@$ptrparamsin),") <> $nomsubsub(",join(",",@$ptrparamsinsubsub),")\n";
              }
              else
              { my $res = test_paramsin2($ptrparamsin,$ptrparamsinsubsub);
                unless($res)
                {
                  print "for subroutine $nomsub, declared in file $file, its call in the subroutine $nomsub2 line $lignesubsub of file $file2";
                  print " is not correct :\n$nomsub(",join(",",@$ptrparamsin),") <> $nomsubsub(",join(",",@$ptrparamsinsubsub),")\n";
                }
              }
            }
          }
        }
      }
    }
  }
}
print "\n";
exit;


sub test_paramsin2
{ my ($params,$ptypes) = @_;
  my $ok=1;
  foreach my $type (@$ptypes)
  { my $found=0;
    foreach my $param (@$params)
    { $found++ if ($type eq $param);
    }
    $ok=0 unless ($found);
  }
  return $ok;
}

sub test_paramsin
{ my ($params,$ptypes,$nomsub,$file) = @_;
  foreach my $type (@$ptypes)
  { my $found=0;
    foreach my $param (@$params)
    { $found++ if ($type eq $param);
    }
    unless ($found || $nomsub eq 'specfem3D' || $nomsub eq 'meshfem3D')
    { print "---\nstructure $type not given as parameters in subroutine $nomsub of file $file\n---\n",join(",",@$params),"\n",join(",",@$ptypes),"\n";
    }
  }
  return;
}

sub analyse_sub
{ my ($fichier,$nom,$pline_num,@paramsin) = @_;
  my $ptypes = [];
  my $pcursub = [$nom,\@paramsin,$ptypes];
  my %hnom;
  while (my $line=<DF>)
  { $$pline_num++;
    if ($line =~ m/^\s*type\s*\([ a-zA-Z_0-9-]+\)\s+(.*)$/)
    { my $type = $1;
      chomp($type);
      $type=trim($type);
      push(@$ptypes,$type);
      $htypes{$type}++;
    }
    # appel sur une ligne
    elsif (($line =~ m/^(?:[^!]*\s+|\s*)call\s+([a-zA-Z_0-9-]*)\s*(?:\(([^!]*)\))?\s*$/) && ($line !~ m/MPI_[A-Z]+/))
    { my $nomsub = $1;
      my $numligne = $$pline_num;
#     unless (exists($hnom{$nomsub}))
#     { $hnom{$nomsub}++;
        my $subparams = [];
        if (defined($2))
        { push(@$subparams,map(trim($_),split(/,/,$2)));
          filtre($subparams,$ptypes);
        }
        my @sub_called;
        push(@sub_called, $nomsub,$numligne, $subparams);
        push (@$pcursub,\@sub_called);
#     }
    }
    # appel sur plusieurs lignes (&)
    elsif (($line =~ m/^(?:[^!]*\s+|\s*)call\s+([a-zA-Z_0-9-]*)\((.*)&\s*/) && ($line !~ m/MPI_[A-Z]+/))
    { my $nomsub = $1;
      my $numligne = $$pline_num;
#     unless (exists($hnom{$nomsub}))
#                 {       $hnom{$nomsub}++;
        my $subparams = [];
        push(@$subparams,map(trim($_),split(/,/,$2)));
                          pop(@$subparams) if ($subparams->[-1] =~ m/^\s*$/);
                          while(<DF>)
                          {       $$pline_num++;
                                  if ($_ =~ m/^\s*(.*)&\s*$/)
                                  {       push(@$subparams,map(trim($_),split(/,/,$1)));
                                          pop(@$subparams) if ($subparams->[-1] =~ m/^\s*$/);
                                  }
                                  elsif ($_ =~ m/^\s*(.*)\)\s*$/)
                                  {       push(@$subparams,map(trim($_),split(/,/,$1)));
                                          pop(@$subparams) if ($subparams->[-1] =~ m/^\s*$/);
                                          last;
                                  }
                                  else
                                  {       print "probleme :  appel non conforme dans $fichier, sub ",$pcursub->[0]," l $$pline_num\n";exit;
                                  }
                          }
                          filtre($subparams,$ptypes);
        my @sub_called;
                          push(@sub_called, $nomsub,$numligne, $subparams);
                          push (@$pcursub,\@sub_called);
#     }
    }
    elsif ($line =~ m/^\s*end\s*subroutine\s+$nom.*$/)
    { last;
    }
    else {}
  }
  # retirer de paramsin les params n'etant pas dans type
  filtre(\@paramsin,$ptypes);
  return $pcursub;
}

sub filtre
{ my ($pparams,$ptypes) = @_;
  my @params_filtered;
  foreach $one_par (@$pparams)
  { my $found=0;
    foreach $one_type (@$ptypes)
    { if ($one_par eq $one_type)
      { $found=1;
        last;
      }
    }
    push(@params_filtered, $one_par) if ($found);
  }
  @$pparams = @params_filtered;
  return;
}

sub trim
{ my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}

