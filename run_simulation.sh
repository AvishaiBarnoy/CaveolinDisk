#!/bin/bash

# script to run optimization.
# 1. run generate_geom.py to create folders and initial geometry files.
# 2. run this script that uses optimize.py

for i in dn*
do
  cd $i
  cp ../optimize.py .
  for j in geom_*
  do
    cd $j
    cp ../optimize.py .
    python optimize.py geom_*.log
    rm -f optimize.py
    cd ..
  done
  cd ..
done
