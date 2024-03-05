# Caveolin-1 in Caveolae simulation
## About
2D simulation of caveolin-1 on caveolae surface

## Description
caveolin is modelled as a disk and are connected to each other by membrane.
Differnet combination of connectiong between neighborig caveolin are computed, and the geometry is minimized with respect to an ideal angle between edges.

## How to use
### Running simulation
In the `generate_geom.py` file edit the following parameters:
- `n_disks` number of proteins in calculations.
- `L` half-distance between proteins in simulation.
- `r_disk` radius of proteins.
1. Run `generate_geom.py` to generate geometries and place them in correct folders
2. Run `run_simulation.sh` to go through folders and optimize each geometry
### Analysis 
To be added in the future

## Development
In the future I will add an input file and rearrange some of the functions.

## Thanks
Thanks for Itai Knaan-Harpaz for the algorithm for the initial geometry
