#!/bin/bash
echo "Prepare interpreter"
if [ $RUNNER_OS == 'Linux' ]; then
  echo "Installing deps for linux"
  echo "... none for now"
  # NOTE: leidenalg installs igraph (in a somewhat messy way for now
  #sudo add-apt-repository -y ppa:igraph/ppa
  #sudo apt-get -qq update
  #sudo apt-get install igraph
elif [ $RUNNER_OS == 'macOS' ]; then
  echo "Find out OSX version"
  osx_version=$(sw_vers -productVersion)
  echo "OSX version: $osx_version"

  echo "Installing deps for OSX"
  # Prepare to exit upon failure
  set -e
  CONDA_URL="https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
  wget -nv "${CONDA_URL}"
  bash Miniconda3-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda

  echo "$PATH"
  export PATH="$HOME/miniconda/bin:$PATH"
  source $HOME/miniconda/bin/activate

  # Make conda environment and activate
  conda create -y -n ci python=$CONDA_PY
  conda activate ci

  # Use pip from conda
  conda install -y pip
  pip --version

else
  echo "OS not recognized: $RUNNER_OS"
  exit 1
fi

# We do not need to actually install anything when deploying because
# it is a pure python package
echo "Installing Python dependencies"
pip install pytest
pip install pytest-cov
pip install coveralls

pip install numpy
pip install scipy
pip install pandas
pip install scikit-learn
pip install loompy
pip install anndata

echo "Install python-igraph. It takes care of installing the igraph C library"
pip install python-igraph

echo "Install leidnalg"
pip install leidenalg

if [ $USE_SCANPY == 'yes' ]; then
 pip install scanpy
fi

# TODO: add scanpy conditionally
