#!/bin/bash

SCRIPT=lowest_dn.py
mkdir analysis

# empty target file
> analysis/e_values.txt

for i in dn*  # {0..7}
do 
  cd $i
  cp ../$SCRIPT .
  pwd
  python $SCRIPT 
  echo ""
  rm -f $SCRIPT 
  cd ..
done
