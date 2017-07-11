#!/bin/bash
set -o errexit -o nounset
PKG_REPO=$PWD
COMMIT="${TRAVIS_COMMIT:-$APPVEYOR_REPO_COMMIT}"
cd ..

echo "Creating repo folder"
mkdir repo
cd repo

## Set up Repo parameters
echo "Git setup."
git init
git config user.name "DeclareDesign Travis"
git config user.email "team@declaredesign.org"
git config push.default simple

## Get drat repo
echo "Upstream drat repo"
git remote add upstream "https://$GH_TOKEN@github.com/DeclareDesign/declaredesign.github.io.git"
git fetch upstream
git checkout master

echo "R script updating repo contents"
Rscript update_repo.R

echo "Add and commit"
git add *
git commit -m "Travis update $PKG_REPO build $COMMIT"
git push

echo "Complete"
