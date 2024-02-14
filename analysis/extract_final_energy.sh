#!/bin/bash
  

echo "combination, geom_energy" > final_geometry_energy.csv
for i in geometries/geom*
do
  line=`head -1 $i/log_energy_*.txt`
  # echo $line
  combination=`echo "$line" |  grep -o '\[.*\]'`
  # echo $combination
  final=`tail -n 1 $i/log_energy_*.txt | awk '{print $2}' | sed 's/e\=//g'`
  echo \"$combination\", $final >> final_geometry_energy.csv  
done
