#!/bin/bash

check_consistency() {
  for i in geometries/geom*
  do
    cp length_consistent.py $i
    cd $i
    python length_consistent.py log_geom*.txt
    rm -f length_consistent.py
    cd ../../
  done
}
check_consistency
