#!/bin/bash

for i in {0..7}
do 
  cd dn0$i
  cp ../lowest_dn.py
  pwd
  python lowest_dn.py
  echo ""
  rm -f lowest_dn.py
  cd ..
done
