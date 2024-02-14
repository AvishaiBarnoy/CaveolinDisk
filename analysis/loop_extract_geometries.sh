#!/bin/bash

#for i in {1..9}
#do 
#  mkdir geometries/geom000$i
#  cp ../results/log_energy_000$i.txt ../results/log_geom_000$i.txt geometries/geom000$i
#done

for i in {10..99}
do
  mkdir geometries/geom00$i
  cp ../results/log_*_00$i.txt geometries/geom00$i
done

for i in {100..999}
do 
  mkdir geometries/geom0$i
  cp ../results/log_*_0$i.txt geometries/geom0$i
done

for i in {1000..4250}
do 
  mkdir geometries/geom$i
  cp ../results/log_*_$i.txt geometries/geom$i
done

for i in geometries/*
do
    #extract_first_last_geomtries.sh
  cp extract_first_last_geometries.sh $i
  cd $i
  ./extract_first_last_geometries.sh log_geom_*.txt
  rm -f extract_first_last_geometries.sh
  cd ../..
done
