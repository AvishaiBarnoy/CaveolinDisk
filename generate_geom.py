from setup import Combinatorics 
from setup import Combination
from setup import group_operations 
import numpy as np
import sys
import os
from tqdm import tqdm

def main(n_disks, r_disk, L):
    
    # Generate combinations
    combinations = Combinatorics(n_disks)
    raw_combinations = combinations.generate_combinations()
    sys.stdout.write(f"Generated {len(raw_combinations)} combinations.\n")
    
    # removed combinations not starting with first element
    raw_combinations = [combination for combination in raw_combinations if combination[0] == 1] 

    # remove combinations using cyclic filter
    filtered_combinations = group_operations.cyclic_filter(raw_combinations, n_disks)
    sys.stdout.write(f"Filtred redundant combinations, now we only have {len(filtered_combinations)} combinations.\n")
    
    # filtered_combinations = filtered_combinations[i_start:i_end]
    for i, combination in tqdm(enumerate(filtered_combinations), total=len(filtered_combinations)):
        combination = Combination(combination, n_disks)
        length_list = combination.map_combination_to_lengths()
        mod_lengths = combination.modify_one_length(disk_radius=7, L=2)
        
        # check if combination is valid
        if combination.is_valid_inequality(L=2) == False:
            continue
        
        points = np.array(combination.calculate_circle_points())

        geom_file = f"geom_{i}.log"
        
        dn = n_disks - len(combination.lengths)
        # print(combination.mod_length)

        # make directories for geometries
        os.makedirs(f'dn{dn}/geom_{i}', exist_ok=True)

        with open(f"dn{dn}/geom_{i}/{geom_file}", "w") as f:
            first_line = f"# combination: {combination.mod_length}\n"  
            f.write(first_line)
            np.savetxt(f, points)

if __name__ == "__main__":
    n_disks = 11
    r_disk = 7
    L = 2
    main(n_disks, r_disk, L)
