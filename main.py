# here should be logic to call different parts
# 1. parse input file
# 2. generate geometries
# 3. adjust analysis scripts
# 4. ask if wants to run minimization
#   4.1 yes - run simulation
#   4.2 no - exit and say goodbye nicely
# 5. auto analyze
# TODO: explain logic of code, what each file does.
# TODO: option for an inputfile generation script that generates all relevant files and scripts 

from parse_inp import read_inputfile
