#!/bin/bash
if [ "$TRAVIS_OS_NAME" == 'osx' ]; then
  export PATH="$HOME/miniconda/bin:$PATH"
  source $HOME/miniconda/bin/activate
  conda activate travis
  PYTHON=$(which python)
  PYTEST="$HOME/miniconda/bin/pytest"
else
  PYTHON=${PYTHON:-python}
  PYTEST=${PYTEST:-"pytest -rxXs --cov=northstar/"}
fi

echo "python: ${PYTHON}"

echo 'Running pytests...'
# LOCAL TESTING:
# PYTHONPATH=$(pwd)/packages:$(pwd):PYTHONPATH pytest -rxXs test

${PYTEST} "test"
