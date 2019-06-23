#!/bin/bash
# only deploy builds for a release_<sematic-version>_RC?? tag to testpypi
if [ -z $TRAVIS_TAG ]; then
  echo 'No TRAVIS_TAG, exit'
  exit 0
fi
TAG1=$(echo $TRAVIS_TAG | cut -f1 -d_)
TAG2=$(echo $TRAVIS_TAG | cut -f2 -d_)
TAG3=$(echo $TRAVIS_TAG | cut -f3 -d_)
if [ -z $TAG2 ]; then
  echo 'No TAG2, exit'
  exit 0;
fi
if [ $TAG1 != 'release' ] || [ "version = $TAG2" != $(cat semiannotate/version.py) ]; then
  echo 'No release tag or wrong version, exit'
  exit 0;
fi
SINGLET_VERSION=$TAG2

# For now deploy only Linux Python 3.7 (it's pure Python)
if [ $TRAVIS_OS_NAME != 'linux' ]; then
  echo 'Not linux, exit'
  exit 0
fi
if [ $TRAVIS_PYTHON_VERSION != '3.7' ]; then
  echo "Not Python 3.7, exit"
  exit 0
fi

# deploy onto pypitest unless you have no RC
if [ -z $TAG3 ]; then
  TWINE_PASSWORD=${TWINE_PASSWORD_PYPI}
  TWINE_REPOSITORY='https://upload.pypi.org/legacy/'
  echo 'Deploying to production pypi'
elif [ ${TAG3:0:2} == 'RC' ]; then
  TWINE_PASSWORD=${TWINE_PASSWORD_TESTPYPI}
  TWINE_REPOSITORY='https://test.pypi.org/legacy/'
  echo 'Deploying to testpypi'
else
  echo "Tag not recognized: $TRAVIS_TAG"
  exit 1
fi
   
echo "Deploying to pip using twine"
echo "TWINE_REPOSITORY=$TWINE_REPOSITORY"
echo "TWINE_USERNAME=$TWINE_USERNAME"
echo "TWINE_PASSWORD=$TWINE_PASSWORD"
pip --version
pip install twine

# Build source
python setup.py sdist --dist-dir dist/

# Upload source
twine upload --repository-url "${TWINE_REPOSITORY}" -u "${TWINE_USERNAME}" -p "${TWINE_PASSWORD}" dist/singlet-${SINGLET_VERSION}.tar.gz

