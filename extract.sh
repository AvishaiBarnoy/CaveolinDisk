#!/bin/bash

SCRIPT=lowest_dn.py
mkdir analysis

# empty target file
> analysis/e_values.txt

for i in dn*
do 
  cd $i
  cp ../$SCRIPT .
  pwd
  python $SCRIPT -e new
  echo ""
  rm -f $SCRIPT 
  cd ..
done
