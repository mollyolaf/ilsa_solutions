#!/bin/bash

# Render the site
quarto render

# Force-refresh the deploy folder
rm -rf _deploy/*
cp -r _site/. _deploy/

# Commit and push
cd _deploy
git add .
git commit -m "Deploy site update: $(date)"
git push origin gh-pages --force
cd ..
