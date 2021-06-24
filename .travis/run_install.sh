#!/bin/bash

# fortran/openMPI compiler
sudo apt-get install -y gfortran libgomp1 openmpi-bin libopenmpi-dev

# python script needs numpy
echo "Python on path: `which python`"
echo "pip on path: $(which pip)"
# for distribution precise
#sudo apt-get install -qq python-numpy #python-scipy
# trusty:
# (Sep2017 update: adding flag --user : see https://github.com/travis-ci/travis-ci/issues/8382)
pip install --user --upgrade pip setuptools wheel
pip install --user --only-binary=numpy numpy
# version info
python --version

# version infos
echo "compiler versions:" ${FC} ${MPIFC} ${CC}
${FC} --version
${MPIFC} --version
${CC} --version

# installs the CUDA toolkit
if [ "$CUDA" == "true" ]; then
  ## distribution precise: from ubuntu 12.04
  #UBUNTU_VERSION=ubuntu1204
  ## distribution trusty: from ubuntu 14.04
  #UBUNTU_VERSION=ubuntu1404
  ## distribution xenial: from ubuntu 16.04
  UBUNTU_VERSION=ubuntu1604

  # CUDA_VERSION - specifies CUDA toolkit version
  ## trusty
  #CUDA_VERSION=6.5-14
  ## xenial
  CUDA_VERSION=9.2.148-1

  echo "Installing CUDA library"
  echo "CUDA version: ${CUDA_VERSION}"
  echo "UBUNTU version: ${UBUNTU_VERSION}"

  # note: travis could stall and time out here
  #       one could try to add: travis_retry sudo dpgk -i ..
  #       https://docs.travis-ci.com/user/common-build-problems/#travis_retry
  #
  # remove old nvidia-cuda packages
  #sudo apt-get remove nvidia-cuda-* ;

  # gets packages
  INSTALLER=cuda-repo-${UBUNTU_VERSION}_${CUDA_VERSION}_amd64.deb
  wget http://developer.download.nvidia.com/compute/cuda/repos/${UBUNTU_VERSION}/x86_64/${INSTALLER}
  sudo dpkg -i ${INSTALLER}

  # ubuntu16.04 version package needs key
  sudo apt-key adv --fetch-keys http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/7fa2af80.pub

  # update
  echo "Updating libraries"
  sudo apt-get update -qq
  dpkg -l | grep cuda
  export CUDA_APT=${CUDA_VERSION:0:3}
  export CUDA_APT=${CUDA_APT/./-}
  echo "CUDA: ${CUDA_APT}"

  # installs packages
  # CUDA_PACKAGES="cuda-drivers cuda-core-${CUDA_APT} cuda-cudart-dev-${CUDA_APT} cuda-cufft-dev-${CUDA_APT}";
  CUDA_PACKAGES="cuda-drivers cuda-compiler-${CUDA_APT} cuda-cudart-dev-${CUDA_APT}"
  echo "Installing ${CUDA_PACKAGES}"
  sudo apt-get install -y --no-install-recommends ${CUDA_PACKAGES}
  sudo apt-get clean
  export CUDA_HOME=/usr/local/cuda-${CUDA_VERSION:0:3}
  export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}
  export PATH=${CUDA_HOME}/bin:${PATH}
  echo ""
  nvcc --version
else
  export CUDA_HOME=""
fi

# storing updated environment parameters for following bash-script
echo "export PATH=${PATH}" > $HOME/.tmprc
echo "export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >> $HOME/.tmprc
echo "export CUDA_HOME=${CUDA_HOME}" >> $HOME/.tmprc

# to avoid mpi issues on travis
# see: https://github.com/open-mpi/ompi/issues/1393#issuecomment-187980473
#      https://github.com/open-mpi/ompi/issues/4948
echo "export export OMPI_MCA_btl_vader_single_copy_mechanism=none" >> $HOME/.tmprc
echo "export export OMPI_MCA_btl=^openib" >> $HOME/.tmprc

echo "exports:"
export
echo ""
