#!/bin/bash
set -o errexit -o nounset
PKG_REPO=$PWD
cd ..

addToDrat(){
  mkdir drat
  cd drat

  ## Set up Repo parameters
  git init
  git config user.name "DeclareDesign Travis"
  git config user.email "team@declaredesign.org"
  git config --global push.default simple

  ## Get drat repo
  git remote add upstream "https://$GH_TOKEN@github.com/DeclareDesign/declaredesign.github.io.git"
  git fetch upstream 2>err.txt
  git checkout master

  Rscript -e "path <- ifelse(.Platform\$OS.type == 'windows', file.path('..', '${APPVEYOR_PROJECT_NAME:-$PKG_REPO}'), file.path('..')); \
  for(pkg in dir(path, pattern = ifelse(.Platform\$OS.type == 'windows', '.zip', '.t*z'))) { print(paste('processing', pkg)); \
  drat::insertPackage(file = file.path(path, pkg), \
  repodir = '.', \
  commit='Travis update ${APPVEYOR_PROJECT_NAME:-$PKG_REPO} build ${APPVEYOR_REPO_COMMIT:-$TRAVIS_COMMIT}', branch = 'master') }"

  git push 2>err.txt

}

addToDrat
