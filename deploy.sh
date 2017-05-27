#!/bin/bash
set -o errexit -o nounset
PKG_REPO=$PWD
cd ..

mkdir drat
cd drat

## Set up Repo parameters
git init
git config user.name "DeclareDesign Travis"
git config user.email "team@declaredesign.org"
git config --global push.default simple

## Get drat repo
git remote add upstream "https://$GH_TOKEN@github.com/DeclareDesign/declaredesign.github.io.git"
git fetch upstream
git checkout master

## move file to the right place (/bin/contrib etc)

## do write_PACKAGES which updates the PACKAGES file appropriately

Rscript ../update_repo.R

git add *

git commit -m "Travis update $PKG_REPO build $TRAVIS_COMMIT"

git push
