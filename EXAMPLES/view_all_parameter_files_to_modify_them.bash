#!/bin/bash

# I get the list of files using " ls */Par_fi* */*/Par_fi* */*/*/Par_fi* */*/*/*/Par_fi* | /bin/grep -v "~" | tr '\n' ' ' "
# and I add ../DATA/Par_file manually

#
# Note that if you want to change the name of a parameter or anything else in all the Par_files below in an automatic way
# you can also use a "find" command combined with a "sed" command, here is an example:
#
#     find . -type f -name "*Par_file*" -exec sed -i 's/ATTENUATION_VISCOELASTIC_SOLID/ATTENUATION_VISCOELASTIC/g' {} \;
#
# However this will *not* make the change in ../DATA/Par_file, which you will need to do separately.
#

if [ -z "$EDITOR" ]
then
  EDITOR=vi
fi

$EDITOR ../DATA/Par_file benchmarks/attenuation_benchmark_GJI_2002_versus_normal_modes/Par_file benchmarks/prem_ani_noocean_Q_Cowling_Bolivia_9s/SEMD/DATA/Par_file benchmarks/prem_ani_noocean_Q_Cowling_Vanuatu/SEMD/DATA/Par_file benchmarks/prem_iso_noocean_Q_Cowling_Vanuatu/SEMD/DATA/Par_file benchmarks/undo_attenuation_benchmark_GJI_2016/DATA/Par_file global_PREM_kernels/amplitude/DATA/Par_file global_PREM_kernels/traveltime/DATA/Par_file global_s362ani_shakemovie/DATA/Par_file global_small/DATA/Par_file noise_examples/global_long/DATA/Par_file_NOISE_1_attenuation noise_examples/global_long/DATA/Par_file_NOISE_1_noattenuation noise_examples/global_long/DATA/Par_file_NOISE_2_attenuation noise_examples/global_long/DATA/Par_file_NOISE_2_noattenuation noise_examples/global_long/DATA/Par_file_NOISE_3_noattenuation noise_examples/global_short/DATA/Par_file_NOISE_1_attenuation noise_examples/global_short/DATA/Par_file_NOISE_1_noattenuation noise_examples/global_short/DATA/Par_file_NOISE_2_attenuation noise_examples/global_short/DATA/Par_file_NOISE_2_noattenuation noise_examples/global_short/DATA/Par_file_NOISE_3_noattenuation noise_examples/regional/DATA/Par_file_NOISE_1_attenuation noise_examples/regional/DATA/Par_file_NOISE_1_noattenuation noise_examples/regional/DATA/Par_file_NOISE_2_attenuation noise_examples/regional/DATA/Par_file_NOISE_2_noattenuation noise_examples/regional/DATA/Par_file_NOISE_3_noattenuation noise_examples/test_global/DATA/Par_file_NOISE_1_attenuation noise_examples/test_global/DATA/Par_file_NOISE_1_noattenuation noise_examples/test_global/DATA/Par_file_NOISE_2_attenuation noise_examples/test_global/DATA/Par_file_NOISE_2_noattenuation noise_examples/test_global/DATA/Par_file_NOISE_3_noattenuation noise_examples/test_regional/DATA/Par_file_NOISE_1_attenuation noise_examples/test_regional/DATA/Par_file_NOISE_1_noattenuation noise_examples/test_regional/DATA/Par_file_NOISE_2_attenuation noise_examples/test_regional/DATA/Par_file_NOISE_2_noattenuation noise_examples/test_regional/DATA/Par_file_NOISE_3_noattenuation point_force/DATA/Par_file regional_Greece_noise_small/DATA/Par_file regional_Greece_small/DATA/Par_file regional_MiddleEast/DATA/Par_file small_benchmark_run_to_test_more_complex_Earth/DATA/Par_file small_benchmark_run_to_test_very_simple_Earth/DATA/Par_file
