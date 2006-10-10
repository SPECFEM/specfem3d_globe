#!/usr/bin/perl

use Time::Local;
use Getopt::Std;
use POSIX;

sub Usage{
  print STDERR <<END;

  process_data_new.pl -m CMTFILE -l Start/End -t Ts/Tl -i -p -x Ext -d OutDir
                      -s Sps -P n/p -T Taper_width -c -y t1/v1/... 
                      data_sac_files
  with

       -m -- use CMTSOLUTION file to set event info
       -l Start/End -- start and end of record from o
       -t Ts/Tl -- shortest and longest period
       -i -- transfer record to displacement
       -p -- pick event p and s first arrival into sac headers (t1 and t2)
       -x Ext -- extension of new file name
       -d Dir -- directory to put the output files (default .)
       -s Sps -- sample per second (default 1.0)
       -P n/p -- number of poles and passes for butterworth filter (default 4/2)
       -T taper_width -- taper width to use in 'taper' command of SAC
       -c check sac output on the screen, otherwise print out summary
       -y t0/vel/t0/vel... --- pick times into headers t3/t4/... (rayleigh = 0/3.8; love = 0/4.0)

       data_sac_files --- names of the data sac files to be processed

 Notice:
  1. We require that polezero files are in the same directory as the sac 
     data files which is generally satisfied. All the needed info are 
     taken from sac headers
  2. Origin time is set to PDE + time_shift (given by the CMTSOLUTION)
  3. The displacement outputs after -i option are in the unit of meters
  4. For BH? components, please set -s 20, otherwise interpolation to
     1 sample/second will be performed

END
  exit(1);
}

if (@ARGV == 0) { Usage(); }

if (!getopts('m:l:t:ipx:d:s:P:T:y:c')) {die('Check input arguments\n');}
if ($opt_t) {($tmin, $tmax) = split(/\//,$opt_t);
	     $f1 = 1./$tmax;  $f2=1./$tmin;}
if ($opt_d) {
  $out_dir=$opt_d; if (not -d $out_dir) {mkdir($out_dir,0777);}
}
if ($opt_x) {$ext='.'.$opt_x;}else {$ext="";}
if (!$opt_s) {$dt=1.0;} else {$dt = 1.0/$opt_s;}
if (!$opt_P) {$poles=4;$pass=2;}
else{($poles,$pass)=split(/\//,$opt_P);
     if(not defined $pass or $pass < 1) {$pass=2;}
}
if (!$opt_T) {$opt_T = 0.05;}
if ($opt_l) {($lmin,$lmax) = split(/\//,$opt_l);} else {$lmin = 0; $lmax = 3600;}

$saclst="/opt/seismo-util/source/saclst/saclst";
$phtimes="/opt/seismo-util/bin/phtimes";
if (! -e $saclst) {die(" No $saclst file\n");}
if (! -e $phtimes) {die("No $phtimes file\n");}

$eps=1e-5; $undef=-12345.0; $cundef="-12345";

foreach $file (@ARGV) {
  if (! -f $file) {die("No such file to be processed!!\n"); }
  print "Processing data file $file \n";
  if (not $opt_d) {$outfile = $file.$ext;}
  else { ($filename) = split(" ",`basename $file`);
	 $outfile = "$out_dir/${filename}${ext}";}
  ($filedir) = split(" ",`dirname $file`);
  if ($ext or $opt_d) {system("\\cp -f $file $outfile");}

  # process data
  if ($opt_c) {open(SAC,"|sac2000 ") || die("Can't open sac2000\n");}
  else {open(SAC,"|sac2000 > /dev/null") || die("Can't open sac2000\n");}
  print SAC "echo on\n";
  print SAC "r $outfile\n";

  # add event information from CMTSOLUTION FILE
  if ($opt_m) {
    print "    Adding event info from $opt_m\n";
    ($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tshift,undef,
     $elat,$elon,$edep,undef) = get_cmt($opt_m);
    ($oyear1,$ojday1,$ohr1,$omin1,$osec1,$omsec1)
      =tdiff($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tshift);
    print SAC "ch o gmt $oyear1 $ojday1 $ohr1 $omin1 $osec1 $omsec1\n";
    print SAC "evaluate to tmp 0 - &1,o\n";
    print SAC "ch allt %tmp% iztype io\n";
    print SAC "ch evla $elat evlo $elon evdp $edep \n";
    print SAC "w $outfile\nquit\n";
    close(SAC);
    if ($opt_c) {open(SAC,"|sac2000 ") || die("Can't open sac2000\n");}
    else {open(SAC,"|sac2000 > /dev/null") || die("Can't open sac2000\n");}
    print SAC "echo on\n";
    print SAC "r $outfile\n";}

  if ($opt_l){  # cut record and rtrend and rmean
    print "    Cut the record from o+$lmin to o+$lmax\n";
    (undef,$tmp_o)=split(" ",`$saclst o f $outfile`);
    if (abs($tmp_o-$undef) < $eps) {die("O has not been set to cut record\n");}
    print SAC "setbb begin ( max &1,b ( &1,o + $lmin ) ) \n";
    print SAC "setbb end   ( min &1,e ( &1,o + $lmax ) ) \n";
    print SAC "cut %begin% %end% \n";
    print SAC "r $outfile\n cut off\n  w over \nquit\n";
    close (SAC);
    if ($opt_c) {open(SAC,"|sac2000 ") || die("Can't open sac2000\n");}
    else {open(SAC,"|sac2000 > /dev/null") || die("Can't open sac2000\n");}
    print SAC "echo on\n";
    print SAC "r $outfile\n";}

    print SAC "interp delta $dt\n";

  if ($opt_t) {# filter records
    print "    Filter record at periods $tmin and $tmax\n";
    print SAC "rtrend\n rmean\n taper width $opt_T\n";
    printf SAC ("bp n %4d p $pass cor %12.6f %12.6f\n",$poles,$f1,$f2);
    printf SAC " rtrend\n rmean\n taper width $opt_T\n";}

  if ($opt_i) {  # transfer instrument response
    if (! $opt_t) {die(" Specify bandpass range by -t\n");}
    $f0=$f1*0.8;
    $f3=$f2*1.2;
    (undef,$network,$sta,$comp,$khole)=split(" ",`$saclst knetwk kstnm kcmpnm khole f $outfile`);
    if ($network eq $cundef or $sta eq $cundef or $comp eq $cundef or $khole eq $cundef) {
      die("No network station name or component name or khole defined in $outfile\n");}
    $pzfile=`ls -1 $filedir/SAC_PZs_${network}_${sta}_${comp}_${khole}* | head -1`;
    chomp($pzfile); 
    if (!-f $pzfile) {die(" No pzfile $pzfile \n");}
    print "    Transfer instrument response using pz file $pzfile\n";
    printf SAC ("trans from polezero s %20s to none freq %12.6f%12.6f%12.6f%12.6f \n",
		$pzfile,$f0,$f1,$f2,$f3);}

  if ($opt_p) { # add p and s arrival info
    print "    Add first P and S arrival information\n";
    (undef,$gcarc,$edep)=split(" ",`$saclst gcarc evdp f $outfile`);
    if (abs($gcarc-$undef) < $eps or abs($edep-$undef) < $eps) {die("No gcarc and depth info\n");}
    ($Pph,$Ptime)=split(" ",`$phtimes $edep $gcarc P`);
    ($Sph,$Stime)=split(" ",`$phtimes $edep $gcarc S`);
    if (not $opt_m) {$tshift = 0;}
    print SAC "evaluate to tmp1 $Ptime - $tshift\n";
    print SAC "evaluate to tmp2 $Stime - $tshift\n";
    print SAC "ch t1 %tmp1% t2 %tmp2%\n";
    print SAC "ch kt1 $Pph kt2 $Sph\n";}

  if ($opt_y) { # add arrival times to headers
    @numbers=split(/\//,$opt_y); $npairs = floor(@numbers/2);
    if (@numbers != $npairs*2) {die("Enter -y t0/Vel/t0/Vel/...\n");}
    for ($i=0;$i<$npairs;$i++) {
      $int = $numbers[$i*2]; $slope = $numbers[$i*2+1];
      print "   Add arrival time for waves with group velocity: $slope\n";
      (undef,$dist,$begin,$end)=split(" ",`$saclst dist b e f $outfile`);
      if (abs($dist-$undef) < $eps) {die("Not defined dist\n");}
      $h1 = 3+$i;  $h1 = "t$h1"; $k1 = "k$h1"; $v1 = $int + $dist/$slope;
      (undef,$v1,undef) = sort {$a <=> $b} ($begin, $v1 ,$end);
      printf SAC ("ch $h1 %12.2f\n",$v1);
      print SAC "ch $k1 $h1\n";}}

#  print "   write file $outfile\n";
  print SAC "w $outfile\n";
  print SAC "echo off\nquit\n";
  close(SAC);
}
print "  Done! \n";

#**********************************************************
sub mday2jday {
  my($oyear,$omonth,$oday)=@_;
  $omonth = $omonth-1;  #months range from 0..11
  $time_sec = timegm(3,3,3,$oday,$omonth,$oyear);
  @t = gmtime($time_sec);
  my($jday) = $t[7];
  $jday += 1;
  return ($jday);
}

sub tdiff{
  my ($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tadd)=@_;
  $time = timegm($osec, $omin, $ohr, $oday , $omonth-1, $oyear);
  $time += ($tadd +$omsec/1000); #event_time in machine format
  $msec = sprintf("%03.0f",($time - (floor($time)))*1000);
  $time = floor($time);
#event-time:
  my($sec, $min, $hr, $day , $month, $year,$weekday,$jday) = gmtime($time);
  $month += 1;
  $year += 1900;
  $jday +=1;
  return ($year,$jday,$hr,$min,$sec,$msec);
}

sub get_cmt {
  my ($cmt_file)=@_;
  open(CMT, "$cmt_file") or die("Error opening $cmt_file\n");
  @cmt = <CMT>;
  my($pde,$oyear,$omonth,$oday,$ohr,$omin,$osec1)=split(" ",$cmt[0]);
  my($osec,$omsec)=split(/\./,$osec1); $omsec=$omsec*10;
  my(undef,undef,$tshift)=split(" ",$cmt[2]);
  my(undef,undef,$hdur)=split(" ",$cmt[3]);
  my(undef,$elat)=split(" ",$cmt[4]);
  my(undef,$elon)=split(" ",$cmt[5]);
  my(undef,$edep)=split(" ",$cmt[6]);
  my(undef,$Mrr)=split(" ",$cmt[7]);
  my(undef,$Mtt)=split(" ",$cmt[8]);
  my(undef,$Mpp)=split(" ",$cmt[9]);
  my(undef,$Mrt)=split(" ",$cmt[10]);
  my(undef,$Mrp)=split(" ",$cmt[11]);
  my(undef,$Mtp)=split(" ",$cmt[12]);
  close(CMT);
  return ($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tshift,$hdur,$elat,$elon,$edep,$Mrr,$Mtt,$Mpp,$Mrt,$Mrp,$Mtp);
}
