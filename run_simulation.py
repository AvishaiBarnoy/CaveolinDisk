from setup import Combinatorics 
from setup import Combination
from setup import group_operations 
from optimization import Optimization 
from optimization import Animation 
import numpy as np
import matplotlib.pyplot as plt
import sys
from tqdm import tqdm

def main(n_disks, r_disk, L, plot=False):
    # Generate combinations
    combinations = Combinatorics(n_disks)
    raw_combinations = combinations.generate_combinations()
    sys.stdout.write(f"Generated {len(raw_combinations)} combinations.\n")
    
    # removed combinations not starting with first element
    raw_combinations = [combination for combination in raw_combinations if combination[0] == 1] 
     
    # remove combinations using cyclic filter
    filtered_combinations = group_operations.cyclic_filter(raw_combinations, n_disks)
    sys.stdout.write(f"Filtred redundant combinations, now we only have {len(filtered_combinations)} combinations.\n")
    
    # minimize structures and log energies
    n = 1
    i_start = 1368 
    i_end   = 1368 + 1367
    filtered_combinations = filtered_combinations[i_start:i_end]
    for i, combination in tqdm(enumerate(filtered_combinations), total=len(filtered_combinations)):
        combination = Combination(combination, n_disks)
        length_list = combination.map_combination_to_lengths()
        mod_lengths = combination.modify_one_length(disk_radius=7, L=2)
        
        # check if combination is valid
        if combination.is_valid_inequality(L=2) == False:
            continue
        points = np.array(combination.calculate_circle_points())
        if len(points) == 0:
            continue

        sys.stdout.write(f"Lengths: {mod_lengths}\n")
        sys.stdout.write(f"Minimizing energy for structure {n+i_start}\n")

        F_tot = combination.calc_comb_energy(L=2)
        sys.stdout.write(f"Calculated F_tot for single combination, {F_tot}\n")

        sys.stdout.write("Placed points on circle.\n")
        log_geom = f"dt01/log_geom_{n+i_start}.txt"
        log_energy = f"dt01/log_energy_{n+i_start}.txt"

        k_angle = 2 * (0.5 * 10 * 1/np.tanh(L))
        visualize = False


        # ideal_angle = h/a*de*1/np.sqrt(k*k_t)*1/(1/np.tanh(R/l)+0.5/np.tanh(L/l)) # absolute value
        ideal_phi = 0.2 * 1/(1/np.tanh(7/3.5) + 0.5 / np.tanh(1))  # absolute value
        x = np.pi - 3/2 * ideal_phi
        ideal_angle = x + ideal_phi
        
        # add combination at top of files 
        with open(log_geom, 'a') as f:
            f.write(f"# combination: {mod_lengths}\t energy: {F_tot}\n")
        with open(log_energy, 'a') as f:
            f.write(f"# combination: {mod_lengths}\t energy: {F_tot}\n")
        
        sim = Optimization(vertices=points,
                            ideal_distances=mod_lengths,
                            ideal_angles=[ideal_angle] * len(mod_lengths),
                            k_angle=k_angle, k_edges=25, max_steps=1000, dt=0.01,
                            optimizer="gd",
                            log_energy=log_energy,
                            log_geom=log_geom, n_geom=100
                            )
        sim.log_std()


        if visualize:
            animation = Animation(sim)
            animation.start_animation()
        else:
            while not sim.should_stop():
                sim.update_position()
        sys.stdout.write("\n")
        n += 1

if __name__ == "__main__":
    n_disks = 17
    r_disk = 7
    L = 2
    plotting = False
    main(n_disks, r_disk, L, plotting)
