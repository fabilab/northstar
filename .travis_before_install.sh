#!/bin/bash
echo "Prepare interpreter"
if [ $TRAVIS_OS_NAME == 'linux' ]; then
  echo "Installing deps for linux"
  #sudo add-apt-repository ppa:nschloe/swig-backports -y
  #sudo apt-get -qq update
  #sudo apt-get install -y swig3.0
  sudo apt-get install igraph
elif [ $TRAVIS_OS_NAME == 'osx' ]; then
  echo "Installing deps for OSX"
  if [ $PYTHON_VERSION == "2.7" ]; then
    CONDA_VER='2'
  elif [ $PYTHON_VERSION == "3.7" ]; then
    CONDA_VER='3'
  else
    echo "Miniconda only supports 2.7 and 3.7"
  fi
  curl "https://repo.continuum.io/miniconda/Miniconda${CONDA_VER}-latest-MacOSX-x86_64.sh" -o "miniconda.sh"
  bash "miniconda.sh" -b -p $HOME/miniconda
  echo "$PATH"
  export PATH="$HOME/miniconda/bin:$PATH"
  source $HOME/miniconda/bin/activate
  # Use pip from conda
  conda install -y pip
  pip --version

  conda install -c conda-forge -y igraph

else
  echo "OS not recognized: $TRAVIS_OS_NAME"
  exit 1
fi

echo "Installing dependencies"
## setuptools < 18.0 has issues with Cython as a dependency
#pip install Cython
#if [ $? != 0 ]; then
#    exit 1
#fi

# deps #FIXME: do better
pip install pytest
pip install pytest-cov
pip install coveralls

pip install numpy
pip install scipy
pip install igraph-python
