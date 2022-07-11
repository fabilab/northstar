#!/bin/bash
if [ $TRAVIS_OS_NAME == 'osx' ]; then
  export PATH="$HOME/miniconda/bin:$PATH"
  source $HOME/miniconda/bin/activate
  conda activate ci
fi

# We do not need to actually install anything when deploying because
# it is a pure python package
if [ $TRAVIS_BUILD_STAGE_NAME == "Test" ]; then
  pip install -v '.[atlas-fetcher]'
  if [ $? != 0 ]; then
      exit 1
  fi

elif [ $TRAVIS_BUILD_STAGE_NAME == "Deploy" ]; then
  echo "Nothing to install for deploy stage"

elif [ $TRAVIS_BUILD_STAGE_NAME == "Test_deployed" ]; then
  # only test deployed builds for a release_<sematic-version>_RC?? tag to testpypi
  if [ -z $TRAVIS_TAG ]; then
    echo 'No TRAVIS_TAG, exit'
    exit 1
  fi
  TAG1=$(echo $TRAVIS_TAG | cut -f1 -d_)
  TAG2=$(echo $TRAVIS_TAG | cut -f2 -d_)
  TAG3=$(echo $TRAVIS_TAG | cut -f3 -d_)
  echo "TAG1: $TAG1"
  echo "TAG2: $TAG2"
  echo "TAG3: $TAG3"
  
  if [ -z $TAG2 ]; then
    echo 'No TAG2, exit'
    exit 1;
  fi

 # deploy onto pypitest unless you have no RC
 if [ -z $TAG3 ]; then
   echo "Installing remote version from pypi"
   pip install -v 'northstar[atlas-fetcher]'
 elif [ ${TAG3:0:2} == 'RC' ]; then
   echo "Installing remote version from testpypi"
   pip install -v --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple 'northstar[atlas-fetcher]'
 else
   echo "Tag not recognized: $TRAVIS_TAG"
   exit 1
 fi

fi
