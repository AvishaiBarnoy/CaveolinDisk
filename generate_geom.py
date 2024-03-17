from combinations import Combinatorics
from combinations import Combination
from combinations import group_operations
from combinations import calc_n_disks
import numpy as np
import sys
import os
from tqdm import tqdm
import argparse

def main(r_disk, L, distribute_excess_membrane=0):

    n_disks = calc_n_disks(R=r_disk, L=L)

    sys.stdout.write(f"""
N proteins: {n_disks}
protein radius: {r_disk} nm
half-distance: {L} nm
estimated caveolin radius {(2*r_disk+L)*n_disks / (2 * np.pi)} nm\n""")

    # Generate combinations
    if n_disks >= 12:
        sys.stdout.write("""\n
############
# WARNING! #
############
This misses some combinations, e.g., for n=6 it misses [1,3,5] corresponding to lengths [2,2,2]
The other method is good up to n=11, since it is computationally heavy.\n
\n""")
        combinations = Combinatorics(n_disks)
        raw_combinations = combinations.generate_combinations()
        sys.stdout.write(f"Generated {len(raw_combinations)} combinations.\n")

        # removed combinations not starting with first element
        raw_combinations = [combination for combination in raw_combinations if combination[0] == 1]

        # remove combinations using cyclic filter
        filtered_combinations = group_operations.cyclic_filter(raw_combinations, n_disks)
        sys.stdout.write(f"Filtred redundant combinations, now we only have {len(filtered_combinations)} combinations.\n")
    elif n_disks < 12:
        sys.stdout.write("""############
# WARNING! #
############
This method does not filter mirror combinations [1,2,3] and [3,2,1] and only works up to n=11
""")
        # filter duplicate combinations 
        filtered_combinations = set(group_operations.filter_partitions(n_disks))
        sys.stdout.write(f"Filtred redundant combinations, now we only have {len(filtered_combinations)} combinations.\n")

    # filtered_combinations = filtered_combinations[i_start:i_end]
    for i, combination in tqdm(enumerate(filtered_combinations), total=len(filtered_combinations)):
        if n_disks > 12:
            combination = Combination(combination, n_disks)
            length_list = combination.map_combination_to_lengths()
        if n_disks < 12:
            combination = Combination(lengths=combination, n_disks=n_disks)

        combination.modify_one_length(disk_radius=7, L=L)

        # check if combination is valid
        if combination.is_valid_inequality(L=L) == False:
            continue

        dn = n_disks - len(combination.lengths)

        if distribute_excess_membrane == False:
            pass

        elif distribute_excess_membrane == True:
            membrane_excess = L*2 + L*2/(n_disks - dn) * dn
            combination.mod_length = [membrane_excess if length == L*2 else length for length in combination.mod_length]

        elif distribute_excess_membrane == 2:
            sys.stdout.write("NON UNIFORM DISTIBUTION IS NOT IMPLEMENTED YET")
            sys.stdout.write("currently implemented in the optimization level but not in the initial conditions level")
            sys.stdout.write("exiting program without generating geometries")
            sys.exit(1)

        points = np.array(combination.calculate_circle_points())

        geom_file = f"geom_{i}.log"


        # make directories for geometries

        if dn < 10:
            dn_folder = f"dn0{dn}"
        if 10 <= dn < 100:
            dn_folder = f"dn{dn}"

        os.makedirs(f"{dn_folder}/geom_{i}", exist_ok=True)

        with open(f"{dn_folder}/geom_{i}/{geom_file}", "w") as f:
            first_line = f"# combination: {combination.mod_length}\n"
            f.write(first_line)
            np.savetxt(f, points)
    sys.exit(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="generate_geom",
        description="generates initial geometries, sorts them accoring to ΔN in folders",
        epilog="Not all those who physics are lost"
    )
    parser.add_argument("-c", "--conserve_membrane", action="store_true", help="decideds what to do with excess membrane\
                        generated when generating ΔN (combining adjacent proteins)")
    args = parser.parse_args()

    distribute_membrane = args.conserve_membrane # uniform distribution
    r_disk = 7
    L = 1
    main(r_disk, L, distribute_excess_membrane=distribute_membrane)
