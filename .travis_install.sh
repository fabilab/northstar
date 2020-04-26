#!/bin/bash
if [ $TRAVIS_OS_NAME == 'osx' ]; then
  export PATH="$HOME/miniconda/bin:$PATH"
  source $HOME/miniconda/bin/activate
  conda activate travis
fi

# We do not need to actually install anything when deploying because
# it is a pure python package
if [ $TRAVIS_BUILD_STAGE_NAME != "Deploy" ]; then
  pip install -v '.'
  if [ $? != 0 ]; then
      exit 1
  fi
fi
