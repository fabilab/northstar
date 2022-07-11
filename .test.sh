#!/bin/bash
if [ "$TRAVIS_OS_NAME" == 'osx' ]; then
  export PATH="$HOME/miniconda/bin:$PATH"
  source $HOME/miniconda/bin/activate
  conda activate ci
  PYTHON=$(which python)
  PYTEST=$(which pytest)
else
  PYTHON=${PYTHON:-python}
  # coveralls is misbehaving
  #PYTEST=${PYTEST:-"pytest -rxXs --cov=northstar/"}
  PYTEST=${PYTEST:-"pytest"}
fi

echo "python: ${PYTHON}"

echo 'Running pytests...'
# LOCAL TESTING:
# PYTHONPATH=$(pwd)/packages:$(pwd):PYTHONPATH pytest -rxXs test

${PYTEST} "test"
