# here should be logic to call different parts
# 1. parse input file
# 2. generate geometries - generate_geom.py
# 3. adjust run and analysis scripts - adjust.sh
# 4. ask if wants to run minimization - run_simulation.sh
#   4.1 yes - run simulation
#   4.2 no - exit and say goodbye nicely
# 5. extract files - extract.sh

# TODO: move all cli options to parse_inp
# https://realpython.com/command-line-interfaces-python-argparse/
#       1. maybe usie in cli alternative to inputfile.txt
#           1.1 argparse.SUPPRESS
#           1.2 type
#           1.3 choices
#       2. maybe fromfile_prefix_chars="@" as a format for inputfile.txt
#       3. maybe allow_abbrev=False to prevent confusion
#       4. 

from parse_inp import read_inputfile
import argparse
import sys

parser = argparse.ArgumentParser(
    description="script for input file parsing",
    epilog="Not all those who physics are lost"
)
parser.add_argument("-i", "--inputfile", help="path to input file with simulation instructions", default="inputfile.txt")

args = parser.parse_args()
inp_params = read_inputfile(inputfile=args.inputfile)

sys.stdout.write(f"Final parameters {inp_params}\n")
sys.exit(0)
