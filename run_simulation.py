from setup import DiskEnergyCalculator
from optimization import PolygonSpringSimulation
from optimization import PolygonSpringAnimation
import numpy as np
import matplotlib.pyplot as plt
import sys
from tqdm import tqdm

def main(n_disks, r_disk, L, plot=False):
    # Generate combinations
    raw_combinations = DiskEnergyCalculator.ProteinAnalyzer.generate_combinations(n_disks)
    sys.stdout.write(f"Generated {len(raw_combinations)} combinations.\n")
    combinations = DiskEnergyCalculator.ProteinAnalyzer.cyclic_filter(raw_combinations, n_disks)
    sys.stdout.write(f"Filtred redundant combinations, now we only have {len(combinations)} combinations.\n")

    # Modify lengths
    lengths = DiskEnergyCalculator.ProteinAnalyzer.map_many_combinations(combinations, n_disks)
    modified_lengths = DiskEnergyCalculator.ProteinAnalyzer.modify_many_lengths(lengths)
     
    # Optimize structures and log energies
    n = 1
    with open('log.txt', 'a') as log_file:
        for i, combination in tqdm(enumerate(combinations), total=len(combinations)):
            combination = DiskEnergyCalculator.ProteinConfiguration(combination, n_disks)
            length_list = combination.map_combination_to_lengths()
            mod_lengths = DiskEnergyCalculator.ProteinAnalyzer.modify_one_length(length_list, disk_radius=7, L=2)
            
            if DiskEnergyCalculator.ProteinAnalyzer.is_valid_inequality(mod_lengths, L=2) == False:
                continue
            
            points = np.array(DiskEnergyCalculator.calculate_circle_points(mod_lengths))
            if len(points) == 0:
                continue

            sys.stdout.write(f"Lengths: {mod_lengths}\n")
            sys.stdout.write(f"Minimizing energy for structure {n}\n")

            F_tot = DiskEnergyCalculator.calc_comb_energy(modified_lengths[i], n_disks)
            sys.stdout.write(f"Calculated F_tot for single combination, {F_tot}\n")

            sys.stdout.write("Placed points on circle.\n")
            log_geom = f"results/log_geom_{n}.txt"
            log_energy = f"results/log_energy_{n}.txt"

            k_angle = 2 * (0.5 * 10 * 1/np.tanh(L))
            visualize = False


            # ideal_angle = h/a*de*1/np.sqrt(k*k_t)*1/(1/np.tanh(R/l)+0.5/np.tanh(L/l)) # absolute value
            ideal_phi = 0.2 * 1/(1/np.tanh(7/3.5) + 0.5 / np.tanh(1))  # absolute value
            x = np.pi - 3/2 * ideal_phi
            ideal_angle = x + ideal_phi
             
            sim = PolygonSpringSimulation(vertices=points,
                                          ideal_distances=mod_lengths,
                                          ideal_angles=[ideal_angle] * len(mod_lengths),
                                          k_angle=k_angle, k_edges=25, max_steps=1000, dt=0.01,
                                          optimizer="gd",
                                          log_energy=log_energy, log_geom=log_geom
                                          )
            sim.log_std()


            if visualize:
                animation = PolygonSpringAnimation(sim)
                animation.start_animation()
            else:
                while not sim.should_stop():
                    sim.update_position()
            # log_file.write(f"{combination}, {F_tot}, {last_energy}\n")
            sys.stdout.write("\n")
            n += 1

if __name__ == "__main__":
    # n_disks = int(input("Enter the number of disks: ")
    n_disks = 17
    r_disk = 7
    L = 2
    plotting = False
    main(n_disks, r_disk, L, plotting)
