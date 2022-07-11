#!/bin/bash
# only deploy builds for a release_<sematic-version>_RC?? tag to testpypi
TAG1=$(echo $GITHUB_REF_NAME | cut -f1 -d_)
TAG2=$(echo $GITHUB_REF_NAME | cut -f2 -d_)
TAG3=$(echo $GITHUB_REF_NAME | cut -f3 -d_)
echo "TAG1: $TAG1"
echo "TAG2: $TAG2"
echo "TAG3: $TAG3"

if [ -z $TAG2 ]; then
  echo 'No TAG2, exit'
  exit 0;
fi
if [ $TAG1 != 'release' ] || [ "version = $TAG2" != $(cat northstar/_version.py) ]; then
  echo 'No release tag or wrong version, exit'
  exit 0;
fi
VERSION=$TAG2

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
# Turns out twine needs this
pip install typing-extensions
pip install twine


# Build source
python setup.py sdist --dist-dir dist/

# Upload source
twine upload --repository-url "${TWINE_REPOSITORY}" -u "${TWINE_USERNAME}" -p "${TWINE_PASSWORD}" dist/northstar-${VERSION}.tar.gz

