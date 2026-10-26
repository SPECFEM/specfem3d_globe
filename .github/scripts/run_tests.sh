#!/bin/bash
#
# runs a test example case
#

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

WORKDIR=`pwd`
dir=${TESTDIR}
TESTID=${TESTID:-}
TESTCOV=${TESTCOV:-}

## OpenCL
if [ "${OPENCL}" == "true" ]; then
  GPU_RUNTIME=2
  #OMP_NUM_THREADS=1
  # to avoid issues with kernel compilation optimizations,
  # switch the work-group transformation method
  export POCL_WORK_GROUP_METHOD=cbs
  # debugging
  if [ "${DEBUG}" == "true" ]; then
    export POCL_DEBUG=err,warn
    export POCL_EXTRA_BUILD_FLAGS="-cl-opt-disable"
    export POCL_VECTORIZER_FORCE_VECTOR_WIDTH=1
  fi
fi

# info
echo "work directory: $WORKDIR"
echo `date`
echo
echo "**********************************************************"
echo
echo "test directory: $dir"
echo
echo "**********************************************************"
echo

# bash function for checking seismogram output with reference solutions
my_test(){
  echo "######################################################################################################################"
  echo "testing seismograms"
  ln -s $WORKDIR/utils/scripts/compare_seismogram_correlations.py
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/
  if [[ $? -ne 0 ]]; then exit 1; fi
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/ | grep min/max | cut -d \| -f 3 | awk '{print "correlation:",$1; if ($1 < 0.999 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}'
  if [[ $? -ne 0 ]]; then exit 1; fi
  echo "######################################################################################################################"
}

my_kernel_test(){
  # kernel value test - checks rho/kappa/mu kernel value outputs
  echo "######################################################################################################################"
  echo "testing kernel values"
  file_ref=REF_KERNEL/output_solver.txt
  file_out=output.log        # captures the OUTPUT_FILES/output_solver.txt when running solver since IMAIN was set to standard out
  if [ ! -e $file_ref ]; then echo "Please check if file $file_ref exists..."; ls -alR ./; exit 1; fi
  if [ ! -e $file_out ]; then echo "Please check if file $file_out exists..."; ls -alR ./; exit 1; fi
  # gets reference expected kernel values from REF_KERNEL/ folder
  RHO=`grep -E 'maximum value of rho[[:space:]]+kernel' $file_ref | cut -d = -f 2 | tr -d ' '`
  KAPPA=`grep -E 'maximum value of kappa[[:space:]]+kernel' $file_ref | cut -d = -f 2 | tr -d ' '`
  MU=`grep -E 'maximum value of mu[[:space:]]+kernel' $file_ref | cut -d = -f 2 | tr -d ' '`
  ALPHAV=`grep -E 'maximum value of alphav[[:space:]]+kernel' $file_ref | cut -d = -f 2 | tr -d ' '`
  BETAV=`grep -E 'maximum value of betav[[:space:]]+kernel' $file_ref | cut -d = -f 2 | tr -d ' '`

  # need at least rho & kappa (for acoustic kernels)
  if [ "$RHO" == "" ]; then
    echo "  missing reference kernel values: RHO=$RHO | KAPPA=$KAPPA MU=$MU | ALPHAV=$ALPHAV BETAV=$BETAV"
    echo
    exit 1
  else
    echo "  reference kernel values: RHO=$RHO | KAPPA=$KAPPA MU=$MU | ALPHAV=$ALPHAV BETAV=$BETAV"
  fi
  # compares with test output - using a relative tolerance of 0.001 (1 promille) with respect to expected value
  # final test result
  PASSED=0
  # checks rho kernel value
  if [ "$RHO" != "" ]; then
    VAL=`grep -E 'maximum value of rho[[:space:]]+kernel' $file_out | cut -d = -f 2 | tr -d ' '`
    echo "kernel rho   : $VAL"
    echo "" | awk '{diff=ex-val;diff_abs=(diff >= 0)? diff:-diff;diff_rel=diff_abs/ex;print "  value: expected = "ex" gotten = "val" - difference absolute = "diff_abs" relative = "diff_rel; if (diff_rel>0.001){print "  failed"; exit 1;}else{print "  good"; exit 0;} }' ex=$RHO val=$VAL
    if [[ $? -ne 0 ]]; then PASSED=1; fi
  fi
  # checks kappa kernel value
  if [ "$KAPPA" != "" ]; then
    VAL=`grep -E 'maximum value of kappa[[:space:]]+kernel' $file_out | cut -d = -f 2 | tr -d ' '`
    echo "kernel kappa : $VAL"
    echo "" | awk '{diff=ex-val;diff_abs=(diff >= 0)? diff:-diff;diff_rel=diff_abs/ex;print "  value: expected = "ex" gotten = "val" - difference absolute = "diff_abs" relative = "diff_rel; if (diff_rel>0.001){print "  failed"; exit 1;}else{print "  good"; exit 0;} }' ex=$KAPPA val=$VAL
    if [[ $? -ne 0 ]]; then PASSED=1; fi
  fi
  # checks mu kernel value
  if [ "$MU" != "" ]; then
    VAL=`grep -E 'maximum value of mu[[:space:]]+kernel' $file_out | cut -d = -f 2 | tr -d ' '`
    echo "kernel mu    : $VAL"
    echo "" | awk '{diff=ex-val;diff_abs=(diff >= 0)? diff:-diff;diff_rel=diff_abs/ex;print "  value: expected = "ex" gotten = "val" - difference absolute = "diff_abs" relative = "diff_rel; if (diff_rel>0.001){print "  failed"; exit 1;}else{print "  good"; exit 0;} }' ex=$MU val=$VAL
    if [[ $? -ne 0 ]]; then PASSED=1; fi
  fi
  # checks alphav kernel value (if anisotropic kernels)
  if [ "$ALPHAV" != "" ]; then
    VAL=`grep -E 'maximum value of alphav[[:space:]]+kernel' $file_out | cut -d = -f 2 | tr -d ' '`
    echo "kernel alphav: $VAL"
    echo "" | awk '{diff=ex-val;diff_abs=(diff >= 0)? diff:-diff;diff_rel=diff_abs/ex;print "  value: expected = "ex" gotten = "val" - difference absolute = "diff_abs" relative = "diff_rel; if (diff_rel>0.001){print "  failed"; exit 1;}else{print "  good"; exit 0;} }' ex=$ALPHAV val=$VAL
    if [[ $? -ne 0 ]]; then PASSED=1; fi
  fi
  # checks betav kernel value (if anisotropic kernels)
  if [ "$BETAV" != "" ]; then
    VAL=`grep -E 'maximum value of betav[[:space:]]+kernel' $file_out | cut -d = -f 2 | tr -d ' '`
    echo "kernel betav : $VAL"
    echo "" | awk '{diff=ex-val;diff_abs=(diff >= 0)? diff:-diff;diff_rel=diff_abs/ex;print "  value: expected = "ex" gotten = "val" - difference absolute = "diff_abs" relative = "diff_rel; if (diff_rel>0.001){print "  failed"; exit 1;}else{print "  good"; exit 0;} }' ex=$BETAV val=$VAL
    if [[ $? -ne 0 ]]; then PASSED=1; fi
  fi
  # overall pass
  if [[ $PASSED -ne 0 ]]; then
    echo "testing kernel values: failed"; exit 1;
  else
    echo "testing kernel values: all good"
  fi
  echo "######################################################################################################################"
}
# test example
cd $dir

# default setup
if [ ! "${RUN_KERNEL}" == "true" ]; then
  # limit number of time steps
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.5:" DATA/Par_file
  # shortens output interval to avoid timeouts
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file
fi

# specific example setups
if [ "${TESTDIR}" == "EXAMPLES/global_small" ]; then
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.1:" DATA/Par_file
fi
if [ "${TESTDIR}" == "EXAMPLES/regional_Berkeley" ]; then
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 4.0:" DATA/Par_file  # needs increase due to source time function
fi
if [ "${TESTDIR}" == "EXAMPLES/regional_Greece_noise_small" ]; then
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.1:" DATA/Par_file
  sed -i "s:2999:199:g" run_this_example.kernel.sh
  # uses kernel script by default
  cp -v run_this_example.kernel.sh run_this_example.sh
fi

# debug
if [ "${DEBUG}" == "true" ]; then
  # limit for debugging
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
fi

# full gravity
if [ "${FULL_GRAVITY}" == "true" ]; then
  # turns on full gravity
  sed -i "s:^FULL_GRAVITY .*:FULL_GRAVITY = .true.:" DATA/Par_file
  # PETSc
  if [ "${PETSC}" == "true" ]; then
    # switch to PETSc Poisson solver
    sed -i "s:^POISSON_SOLVER .*:POISSON_SOLVER = 1:" DATA/Par_file
  fi
  # set NSTEP for short checks only
  echo "NSTEP = 2" >> DATA/Par_file
fi

## HDF5 - i/o example
if [ "${HDF5}" == "true" ]; then
  echo
  echo "test run w/ HDF5"
  echo
  # turns on HDF5
  echo "turning on HDF5"
  sed -i "s:^HDF5_ENABLED .*:HDF5_ENABLED    = .true.:" DATA/Par_file
  #sed -i "s:^HDF5_FOR_MOVIES .*:HDF5_FOR_MOVIES    = .true.:" DATA/Par_file
  #sed -i "s:^HDF5_IO_NODES .*:HDF5_IO_NODES    = 1:" DATA/Par_file
  # replaces run script
  #cp -v run_this_example_HDF5_IO_server.sh run_this_example.sh
fi

## adios
if [ "${ADIOS2}" == "true" ]; then
  # turns on ADIOS
  echo "turning on ADIOS"
  sed -i "s:^ADIOS_ENABLED .*:ADIOS_ENABLED = .true.:" DATA/Par_file
fi

## GPU
if [ "${GPU}" == "true" ]; then
  # turns on GPU
  echo "turning on GPU"
  sed -i "s:^GPU_MODE .*:GPU_MODE    = .true.:" DATA/Par_file
  sed -i "s:^GPU_PLATFORM .*:GPU_PLATFORM    = *:" DATA/Par_file
  sed -i "s:^GPU_DEVICE .*:GPU_DEVICE    = *:" DATA/Par_file
  # GPU runtime: 0=compile-time, 1=CUDA, 2=OpenCL, 3=HIP
  if [ "${GPU_RUNTIME}" != "" ]; then
    sed -i "s:^GPU_RUNTIME .*:GPU_RUNTIME    = ${GPU_RUNTIME}:" DATA/Par_file
  fi
fi

# coverage runs use short steps
if [ "$TESTCOV" == "true" ]; then
  # limit for debugging
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  if [ "${TESTDIR}" == "EXAMPLES/regional_Greece_noise_small" ]; then
    # will have a number of time steps = 9
    # add a line after generating S_square file to delete lines > 9 (for running dummy simulation)
    sed -i "/run_generate_S_squared/a sed -i '10,$ d' NOISE_TOMOGRAPHY/S_squared" run_this_example.sh
  else
    # set NSTEP for short checks only
    echo "NSTEP = 2" >> DATA/Par_file
  fi
fi

# save Par_file state
if [ -e DATA/Par_file ]; then
  cp -v DATA/Par_file DATA/Par_file.bak
fi

# runs simulation
if [ "${RUN_KERNEL}" == "true" ]; then
  # use kernel script
  ./run_this_example_kernel.sh | tee output.log
else
  # default script
  ./run_this_example.sh
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# simulation done
echo
echo "simulation done: `pwd`"
echo `date`
echo

# seismogram comparison
RUN_COMPARE=true
# turn off for non-default runs
if [ "${TESTCOV}" == "true" ]; then RUN_COMPARE=false; fi
if [ "${DEBUG}" == "true" ]; then RUN_COMPARE=false; fi
if [ "${RUN_KERNEL}" == "true" ]; then RUN_COMPARE=false; fi
if [ "${FULL_GRAVITY}" == "true" ]; then RUN_COMPARE=false; fi
if [ "${TESTDIR}" == "EXAMPLES/regional_Greece_noise_small" ]; then RUN_COMPARE=false; fi

if [ "${RUN_COMPARE}" == "true" ]; then
  my_test
else
  # no comparisons
  :     # do nothing
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# kernel test
if [ "${RUN_KERNEL}" == "true" ]; then
  # check kernel values
  my_kernel_test
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
  # clean up
  rm -rf OUTPUT_FILES/ SEM/ output.log

  # re-run kernel test w/ UNDO_ATT
  UNDO_ATT=`grep ^UNDO_ATTENUATION DATA/Par_file | cut -d = -f 2 | tr -d ' '`
  if [[ ${UNDO_ATT} == *"false"* ]]; then
    echo
    echo "*****************************************"
    echo "run kernel w/ UNDO_ATTENUATION"
    echo "*****************************************"
    echo

    # turns on UNDO_ATTENUATION
    echo "turning on UNDO_ATTENUATION"
    sed -i "s:^UNDO_ATTENUATION .*:UNDO_ATTENUATION = .true.:" DATA/Par_file

    # use kernel script
    ./run_this_example_kernel.sh | tee output.log
    # checks exit code
    if [[ $? -ne 0 ]]; then exit 1; fi
    # kernel test
    my_kernel_test
    # checks exit code
    if [[ $? -ne 0 ]]; then exit 1; fi
  fi
fi

# restore original Par_file
if [ -e DATA/Par_file.bak ]; then
  cp -v DATA/Par_file.bak DATA/Par_file
fi

# cleanup
rm -rf OUTPUT_FILES*
if [ -e DATABASES_MPI ]; then rm -rf DATABASES_MPI*; fi
if [ -e SEM ]; then rm -rf SEM/; fi

echo
echo "all good"
echo `date`
echo
