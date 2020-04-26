#!/bin/bash
echo "Prepare interpreter"
if [ $TRAVIS_OS_NAME == 'linux' ]; then
  echo "Installing deps for linux"
  # NOTE: leidenalg installs igraph (in a somewhat messy way for now)
  #sudo add-apt-repository -y ppa:igraph/ppa
  #sudo apt-get -qq update
  #sudo apt-get install igraph
elif [ $TRAVIS_OS_NAME == 'osx' ]; then
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
  conda create -y -n travis python=$CONDA_PY
  conda activate travis

  # Use pip from conda
  conda install -y pip
  pip --version

else
  echo "OS not recognized: $TRAVIS_OS_NAME"
  exit 1
fi

# We do not need to actually install anything when deploying because
# it is a pure python package
if [ $TRAVIS_BUILD_STAGE_NAME != "Deploy" ]; then
  echo "Installing Python dependencies"
  pip install pytest
  pip install pytest-cov
  pip install coveralls
  
  pip install numpy
  pip install scipy
  pip install pandas
  pip install scikit-learn
  pip install loompy
  
  echo "Install python-igraph. It takes care of installing the igraph C library"
  pip install python-igraph
  
  echo "Install leidnalg"
  pip install leidenalg
  #echo "Install development version of leidenalg"
  #git clone --single-branch --branch develop https://github.com/vtraag/leidenalg.git
  #cd leidenalg
  #python setup.py install
fi
