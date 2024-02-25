#!/bin/bash


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
