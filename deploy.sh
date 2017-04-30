#!/bin/bash
set -o errexit -o nounset
PKG_REPO=$PWD
cd ..

addToDrat(){
  mkdir drat; cd drat

  ## Set up Repo parameters
  git init
  git config user.name "DeclareDesign Travis"
  git config user.email "team@declaredesign.org"
  git config --global push.default simple

  ## Get drat repo
  git remote add upstream "https://$GH_TOKEN@github.com/DeclareDesign/declaredesign.github.io.git"
  git fetch upstream 2>err.txt
  git checkout master

  Rscript -e "for(pkg in dir('..', pattern = '.t*z')) { drat::insertPackage(paste0('../', pkg), \
  repodir = '.', \
  commit='Travis update $PKG_REPO: build $TRAVIS_BUILD_NUMBER', branch = 'master') }"

  git push 2>err.txt

}

addToDrat
