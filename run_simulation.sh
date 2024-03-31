#!/bin/bash

# script to run optimization.
# 1. run generate_geom.py to create folders and initial geometry files.
# 2. run this script that uses optimize.py



CAVEOLIN=`pwd`

# default values for N and L, just to initialize them
N=
L=

# Function to display script usage
usage() {
  echo "Usage: $(basename "$0") [-N value] [-L value] [-h]"
  echo "  -N    Specify value for N, number of proteins"
  echo "  -L    Specify value for L, initial membrane between proteins"
  echo "  -h    Display this help message"
  exit 1
}

# Parse command line options
while getopts "N:L:h" opt; do
  case ${opt} in
    N )
      N=$OPTARG
      ;;
    L )
      L=$OPTARG
      ;;
    h )
      usage
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      usage
      ;;
  esac
done
shift 1
# shift $((OPTIND -1))
# echo $((OPTIND -1))

# Check if no arguments are provided, then display usage
if [ "$#" -eq 0 ]; then
  usage
fi

# if [ "$L" -eq 0 ]; then
if [ -z "$L" ]; then
  echo "Error: L not provided. Please provide L using the -L flag."
  exit 1
fi

# make folder for results
mkdir -p results

# make folder inside results
mkdir -p results/N"$N"/L"$L"
cp -r scripts/* results/N"$N"/L"L"

cd results/N"$N"

# adjust scripts to use correct L value
./adjust_L.sh $L

# generates initial geometries
if [ -z "$N" ]; then
  # if N not provided auto calculate N
  python generate_geom.py -c
else
  # if N is provided pass intp generate_geom
  python generate_geom.py -c -N $N -L $L
fi

# minimizes geometries, adjust this loop to add more/different minimization steps 
for i in dn*/geom_*/*log
do
  OUTPUT_FILE=$(dirname "$i")/geom_opt.txt
  python optimize.py -e new -s -i "$i" -o "$OUTPUT_FILE" -n 25000 -opt cg 
  python optimize.py -e new -s -i "$OUTPUT_FILE" -o "$OUTPUT_FILE" -n 500 -opt bfgs
done

./extract.sh
