#!/bin/bash
if [ $TRAVIS_OS_NAME == 'osx' ]; then
  export PATH="$HOME/miniconda/bin:$PATH"
  source $HOME/miniconda/bin/activate
fi

pip install -v '.'
if [ $? != 0 ]; then
    exit 1
fi
