#!/bin/bash
###################################################

# executable
var=xwrite_profile

# relative location of SPECFEM3D_GLOBE EXAMPLES/ directory for test
EXAMPLES="../../EXAMPLES/global_s362ani_shakemovie"

###################################################

# bash function for checking profile
my_test(){
  echo "testing profile:";
  # remove trailing comments like ... # moho
  sed 's/ #.*$//g' REF_DATA/CARDS_th0042_ph0338 > tmp1.dat
  sed 's/ #.*$//g' OUTPUT_FILES/CARDS_th0042_ph0338 > tmp2.dat
  # remove comment lines starting with # ..
  sed '/^[[:blank:]]*#.*$/d;s/#.//' tmp1.dat > tmpA.dat
  sed '/^[[:blank:]]*#.*$/d;s/#.//' tmp2.dat > tmpB.dat
  # joins file lines
  paste tmpA.dat tmpB.dat > tmp.dat;
  # debug
  #head -n 20 tmp.dat
  # compare radius
  # more sensitive test (works only w/ gfortran tests)
  #awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($1 - $10)**2;print $0," - val = ",val;}}END{print "L2 radius = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  # less sensitive test
  # works w/ intel ifort: two of the ouputted radius locations are shifted by 1km (val > 2.01 check) due to a different numerical precision
  #                       all other outputted values are still the same as w/ gfortran
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($1 - $10)**2;}}END{print "L2 radius = ",val;if(val>2.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then echo "error model: $model "; echo "comparison failed, please check..."; exit 1; fi
  # compare rho
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($2 - $11)**2;}}END{print "L2 rho = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then echo "error model: $model "; echo "comparison failed, please check..."; exit 1; fi
  # compare vpv
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($3 - $12)**2;}}END{print "L2 vpv = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then echo "error model: $model "; echo "comparison failed, please check..."; exit 1; fi
  # compare vsv
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($4 - $13)**2;}}END{print "L2 vsv = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then echo "error model: $model "; echo "comparison failed, please check..."; exit 1; fi
  # compare Qkappa
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($5 - $14)**2;}}END{print "L2 Qkappa = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then echo "error model: $model "; echo "comparison failed, please check..."; exit 1; fi
  # compare Qmu
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($6 - $15)**2;}}END{print "L2 Qmu = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then echo "error model: $model "; echo "comparison failed, please check..."; exit 1; fi
  # compare vph
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($7 - $16)**2;}}END{print "L2 vph = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then echo "error model: $model "; echo "comparison failed, please check..."; exit 1; fi
  # compare vsh
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($8 - $17)**2;}}END{print "L2 vsh = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then echo "error model: $model "; echo "comparison failed, please check..."; exit 1; fi
  # compare eta
  awk 'BEGIN{val=0;}{if(index($0,"#") == 0){val+=($9 - $18)**2;}}END{print "L2 eta = ",val;if(val>0.01){print "failed",val;exit 1;}else{print "good";exit 0;}}' tmp.dat;
  if [[ $? -ne 0 ]]; then echo "error model: $model "; echo "comparison failed, please check..."; exit 1; fi
}

testdir=`pwd`
me=`basename "$0"`

#checks if ROOT valid
if [ -z "${ROOT}" ]; then export ROOT=../../ ; fi

# title
echo >> $testdir/results.log
echo "$me in: $testdir" >> $testdir/results.log
echo >> $testdir/results.log

#cleanup output
rm -rf ./DATABASES_MPI ./OUTPUT_FILES ./tmp.dat
mkdir -p OUTPUT_FILES DATABASES_MPI

# setup
cp -v $EXAMPLES/DATA/Par_file ./DATA/ >> $testdir/results.log

# links should be done by running configure
#cd ./DATA
#ln -s ../../../DATA/crust1.0
#ln -s ../../../DATA/crust2.0
#ln -s ../../../DATA/crustmap
#ln -s ../../../DATA/epcrust
#ln -s ../../../DATA/eucrust-07
#ln -s ../../../DATA/full_sphericalharmonic_model
#ln -s ../../../DATA/s20rts
#ln -s ../../../DATA/s40rts
#ln -s ../../../DATA/s362ani
#ln -s ../../../DATA/sglobe
#ln -s ../../../DATA/QRFSI12
#ln -s ../../../DATA/moon
#ln -s ../../../DATA/topo_bathy
#cd ../


make $var >> $testdir/results.log 2>&1
# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "setup failed, please check..." >> $testdir/results.log
  exit 1
fi

# test models
tests=( 1d_1066a 1d_ak135f_no_mud 1d_ref 1d_jp3d 1d_sea99 1d_ccrem \
        1d_isotropic_prem_onecrust transversely_isotropic_prem_plus_3D_crust_2.0 \
        s362ani s362ani_3DQ \
        s362iso_crustmaps s362wmani_EPcrust s29ea_EUcrust s29ea_crust1.0 \
        s20rts s40rts_crust1.0 sgloberani_aniso full_sh \
        ishii \
        1d_sohl 1d_case65TAY_3D_crust vpremoon )

# single profile tests
for model in ${tests[@]};
do
  echo >> $testdir/results.log
  echo "testing model: $model" >> $testdir/results.log
  echo >> $testdir/results.log

  # sets model
  sed -i "s:^MODEL .*:MODEL   = $model:" DATA/Par_file
  sed -i "s:^ELLIPTICITY .*:ELLIPTICITY   = .false.:" DATA/Par_file    # no ellipticity profiles

  # turns off OCEANS for Mars & Moon models
  if [ "$model" == "1d_sohl" ] || [ "$model" == "1d_case65TAY_3D_crust" ] || [ "$model" == "vpremoon" ]; then
    sed -i "s:^OCEANS .*:OCEANS   = .false.:" DATA/Par_file
  else
    sed -i "s:^OCEANS .*:OCEANS   = .true.:" DATA/Par_file
  fi

  # cleans output
  rm -rf OUTPUT_FILES/*

  # MODEL s362ani profile
  mpirun -np 1 ./bin/xwrite_profile 1 1 47.5 -22.0 >> $testdir/results.log 2>&1

  # checks exit code
  if [[ $? -ne 0 ]]; then
    echo >> $testdir/results.log
    echo "running xwrite_profile failed, please check..." >> $testdir/results.log
    exit 1
  fi

  # sets reference solution
  rm -f REF_DATA/CARDS_th0042_ph0338
  cp -v REF_DATA/CARDS_th0042_ph0338_${model} REF_DATA/CARDS_th0042_ph0338  >> $testdir/results.log
  if [[ $? -ne 0 ]]; then
    echo >> $testdir/results.log
    echo "setup failed, please check..." >> $testdir/results.log
    exit 1
  fi

  # debug
  #echo ""
  #head -n 20 OUTPUT_FILES/CARDS_th0042_ph0338
  #echo ""
  #tail -n 20 OUTPUT_FILES/CARDS_th0042_ph0338
  #echo ""

  # test seismograms
  my_test >> $testdir/results.log

  # checks exit code
  if [[ $? -ne 0 ]]; then
    echo >> $testdir/results.log
    echo "model: $model " >> $testdir/results.log
    echo "run failed, please check..." >> $testdir/results.log
    exit 1
  fi

  # clean card
  rm -f REF_DATA/CARDS_th0042_ph0338 ./tmp*.dat
done

# cleanup
rm -rf ./OUTPUT_FILES ./DATABASES_MPI ./DATA ./tmp*.dat

echo "successful run" >> $testdir/results.log

