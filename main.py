# here should be logic to call different parts
# 1. parse input file
# 2. generate geometries - generate_geom.py
# 3. adjust run and analysis scripts - adjust.sh
# 4. ask if wants to run minimization - run_simulation.sh
#   4.1 yes - run simulation
#   4.2 no - exit and say goodbye nicely
# 5. extract files - extract.sh

from parse_inp import read_inputfile
import argparse
import sys

parser = argparse.ArgumentParser(
    prog="parse_inp",
    description="script for input file parsing",
    epilog="Not all those who physics are lost"
)
parser.add_argument("-i", "--inputfile", help="path to input file with simulation instructions", default="inputfile.txt")

args = parser.parse_args()
inp_params = read_inputfile(inputfile=args.inputfile)

sys.stdout.write(f"Final parameters {inp_params}\n")
sys.exit(0)
