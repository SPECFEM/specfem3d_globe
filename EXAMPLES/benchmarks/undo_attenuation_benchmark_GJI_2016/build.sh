#!/bin/bash

set -e

# Directory of this script; allows to build from somewhere else.
DIR="$(cd $(dirname $0); echo $PWD)"
# Top of specfem3d_globe directory.
BASE="$(cd $DIR/../../..; echo $PWD)"
# Directory : Exact undoing : Partial physical dispersion : Undo attenuation
DATA="\
  exact:.true.:.false.:.false. \
  ppd:.false.:.true.:.false. \
  undo:.false.:.false.:.true. \
"

# Command-line options
while getopts cmh opt; do
  case $opt in
    c)
      do_configure=yes
      ;;
    m)
      do_make=yes
      ;;
    h)
      echo "Usage: $0 [-h] [-c(onfigure-only)] [-m(ake-only)]"
      exit
      ;;
  esac
done
if [ -z "$do_configure" -a -z "$do_make" ]; then
  # No option; do both.
  do_configure=yes
  do_make=yes
fi


# Configure stuff
if [ "x$do_configure" == "xyes" ]; then
for data in $DATA; do
    IFS=: read -r dir exact ppd undo <<< "$data"

  mkdir $dir
  cd $dir

  for sim in 1 3; do
    if [ $sim -eq 1 ]; then
      mkdir step1-forward
      cd step1-forward
      save_forward=".true."
    else
      mkdir step2-adjoint
      cd step2-adjoint
      save_forward=".false."
    fi

    # Configure the build.
    $BASE/configure

    # Set parameters for this study.
    cp -d $DIR/DATA/* DATA/
    sed -i \
      -e "s/^PARTIAL_PHYS_DISPERSION_ONLY.\\+\$/PARTIAL_PHYS_DISPERSION_ONLY = ${ppd}/g" \
      -e "s/^UNDO_ATTENUATION.\\+\$/UNDO_ATTENUATION = ${undo}/" \
      -e "s/^SIMULATION_TYPE.\\+\$/SIMULATION_TYPE = ${sim}/" \
      -e "s/^SAVE_FORWARD.\\+\$/SAVE_FORWARD = ${save_forward}/" \
      DATA/Par_file
    sed -i \
      -e "s/^.\\+EXACT_UNDOING_TO_DISK.\\+\$/  logical, parameter :: EXACT_UNDOING_TO_DISK = ${exact}/" \
      setup/constants.h

    # Create some runtime directories.
    if [ $sim -eq 1 ]; then
      if [ $exact == ".true." ]; then
        mkdir huge_dumps
      fi
      mkdir DATABASES_MPI
    else
      if [ $exact == ".true." ]; then
        ln -s ../step1-forward/huge_dumps
      fi
      ln -s ../step1-forward/DATABASES_MPI
      cd OUTPUT_FILES
      ln -s ../../step1-forward/OUTPUT_FILES/addressing.txt
      cd ..
    fi

    cd ..
  done

  cd ..
done
fi

# Build.
if [ "x$do_make" == "xyes" ]; then
for dir in exact ppd undo; do
  make -j8 -C $dir/step1-forward
  make -j8 -C $dir/step2-adjoint
done
fi
