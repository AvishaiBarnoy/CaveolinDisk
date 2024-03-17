#!/bin/bash

# script to run optimization.
# 1. run generate_geom.py to create folders and initial geometry files.
# 2. run this script that uses optimize.py

# TODO: change everything to run in python
#   or maybe don't copy optimize every time rather run everything from same optimize
#   loop through folder get their path and run optimize and then save output to 
#   same folder.

# TODO: use dirname and basename

for i in dn*
do
  cd $i
  cp ../optimize.py .
  for j in geom_*
  do
    cd $j
    cp ../optimize.py .
    python optimize.py -s -c -i geom_*.log
    rm -f optimize.py
    cd ..
  done
  cd ..
done
