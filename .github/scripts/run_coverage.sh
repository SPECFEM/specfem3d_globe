#!/bin/bash
#
# runs additional coverage examples
#

#set -euo pipefail

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

WORKDIR=$(pwd)
TESTID=${TESTID:-}
TESTCOV=${TESTCOV:-}

if [ "$TESTCOV" != "true" ]; then
  echo "TESTCOV=${TESTCOV} (not coverage), skipping run_coverage.sh"
  exit 0
fi

run_simple() {
  local rel_dir="$1"
  local nstep="$2"
  local mode="$3"

  echo "##################################################################"
  echo "${rel_dir} (simple coverage) NSTEP=${nstep} MODE=${mode}"
  echo

  cd "${WORKDIR}/${rel_dir}"

  # setup
  cp -v DATA/Par_file DATA/Par_file.org
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file
  echo "NSTEP    = ${nstep}" >> DATA/Par_file

  if [ "${mode}" == "GPU" ]; then
    # turns on GPU
    echo "turning on GPU"
    sed -i "s:^GPU_MODE .*:GPU_MODE = .true.:" DATA/Par_file
  elif [ "${mode}" == "FULL_GRAVITY" ]; then
    # turns on full gravity
    sed -i "s:^FULL_GRAVITY .*:FULL_GRAVITY = .true.:" DATA/Par_file
  fi

  # run
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi

  # cleanup
  mv -v DATA/Par_file.org DATA/Par_file
  rm -rf OUTPUT_FILES*
  if [ -e DATABASES_MPI ]; then rm -rf DATABASES_MPI*; fi

  cd "$WORKDIR"
}

run_kernel() {
  local rel_dir="$1"
  local nstep="$2"
  local mode="$3"

  echo "##################################################################"
  echo "${rel_dir} (kernel coverage) NSTEP=${nstep} MODE=${mode}"
  echo

  cd "${WORKDIR}/${rel_dir}"

  # setup
  cp -v DATA/Par_file DATA/Par_file.org
  sed -i "s:^RECORD_LENGTH_IN_MINUTES .*:RECORD_LENGTH_IN_MINUTES = 0.0:" DATA/Par_file
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file
  echo "NSTEP    = ${nstep}" >> DATA/Par_file

  if [ "${mode}" == "GPU" ]; then
    # turns on GPU
    echo "turning on GPU"
    sed -i "s:^GPU_MODE .*:GPU_MODE = .true.:" DATA/Par_file
  fi

  if [ "${rel_dir}" == "EXAMPLES/regional_EMC_model" ]; then
    sed -i "s:^t_start.*:t_start=-4.5:" create_adjoint_sources.sh
    sed -i "s:^t_end.*:t_end=-4.0:" create_adjoint_sources.sh
  fi

  # run
  ./run_this_example_kernel.sh | tee output.log
  if [[ $? -ne 0 ]]; then exit 1; fi

  # cleanup
  mv -v DATA/Par_file.org DATA/Par_file
  rm -rf OUTPUT_FILES*
  if [ -e DATABASES_MPI ]; then rm -rf DATABASES_MPI*; fi
  if [ -e SEM ]; then rm -rf SEM/; fi

  cd "$WORKDIR"
}


echo
echo "coverage run: TESTID=${TESTID}"
echo "work directory: ${WORKDIR}"
echo

# additional example tests (after base to avoid repeating code setup/configuration/compilation)
case "$TESTID" in
  0) # default
    run_simple "EXAMPLES/point_force" 2
    run_simple "EXAMPLES/regular_kernel" 2
    run_simple "EXAMPLES/regional_Greece_small_LDDRK" 2
    ;;
  1) # vectorization
    run_simple "EXAMPLES/regional_s40rts" 2
    run_simple "EXAMPLES/regional_Berkeley" 2
    run_simple "EXAMPLES/mars_regional" 2
    run_simple "EXAMPLES/moon_global" 2
    run_simple "EXAMPLES/regional_sgloberani" 2
    ;;
  2) # vectorization & openMP
    export OMP_NUM_THREADS=2  # limit OpenMP threads
    run_simple "EXAMPLES/global_small" 2
    run_simple "EXAMPLES/global_small" 2 "FULL_GRAVITY"
    ;;
  3) # vectorization w/ ADIOS2
    run_kernel "EXAMPLES/regional_Greece_small_LDDRK" 2
    ;;
  4) # w/ OpenCL
    run_kernel "EXAMPLES/regional_EMC_model" 10 "GPU"
    ;;
  5) # EMC model w/ netCDF & HIP
    run_kernel "EXAMPLES/regional_EMC_model" 10 "GPU"
    ;;
  6) # vectorization w/ NGLL 6
    echo "TESTID=6: no additional coverage examples (base run already executed)"
    ;;
  *)
    echo "TESTID=${TESTID}: no additional coverage examples configured"
    ;;
esac

echo
echo "coverage examples done"
echo "$(date)"
echo
