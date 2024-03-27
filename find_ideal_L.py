from combinations import Combination
import numpy as np
import argparse
from scipy.spatial.distance import euclidean

"""
1. For each L, build geometry
2. Calculate energy for each geometry
3. move to next L
"""

def calculate_new_energy(file_path):
    geometry_energy = 0
    for i in range(len(points)):
        prev_point = points[i - 1]
        current_point = points[i]
        if i == len(points) - 1:
            next_point = points[0]
        else:
            next_point = points[i + 1]

        # get phi_i
        vector1 = prev_point - current_point
        vector2 = next_point - current_point
        dot_product = np.dot(vector1, vector2)
        magnitude_product = np.linalg.norm(vector1) * np.linalg.norm(vector2)
        phi_i= np.arccos(dot_product / magnitude_product)

        # get L_i and R_i
        if i % 2 == 0:
            R_i = euclidean(next_point, current_point) / 2
            L_i = euclidean(prev_point, current_point) / 2
        elif i % 2 != 0:
            R_i = euclidean(prev_point, current_point) / 2
            L_i = euclidean(next_point, current_point) / 2

        geometry_energy += new_energy(L=L_i, phi=phi_i, R=R_i)
    return geometry_energy

def new_energy(L, phi, R, k=0.4e-19, kt=20e-3, h=2, a=0.7, depsilon=4):
    """
    energy calculation for each L and phi
    """
    # constants
    xi = np.sqrt(k/kt) * 1e9 # J / nm
    r = 1/np.tanh(R/xi)
    K = np.sqrt(k*kt) * 1e-9 / 4.11e-21 # kT/nm
    e = h/a * depsilon # kT

    # dependency on L
    l = 1/np.tanh(L/xi)

    # convert phi according to convention
    phi = np.pi - phi

    f1 = 0.5 * K * l * (2*r + l) / (r + l) * phi**2
    f2 = l / (r + l) * e * phi
    f3 = 0.5 / (r + l) * e**2 / K

    F = f1 - f2 - f3

if __name__ == "__main__":
    Li_default = 0.1
    parser = argparse.ArgumentParser(
        description="For each N disks finds optimal half-distance (L) for ΔN=0",
        epilog="Not all those who physics are lost"
    )
    parser.add_argument("-N", "--n_disks", required=True, type=int)
    parser.add_argument("-Li", "--L_initial", default=Li_default)
    parser.add_argument("-Lf", "--L_final", default=Li_default+5)
    args = parser.parse_args()

    # print("N", args.n_disks)

    L = np.arange(args.L_initial, args.L_final, step=0.1)
    energy_lst = np.zeros(L.shape)

    combination = np.array([1]*args.n_disks)
    for i,j in enumerate(L):
        combination = Combination(n_disks=args.n_disks)
        combination.lengths = [1]*args.n_disks # generate initial ΔN=0 
        # combination.map_combination_to_lengths()
        combination.modify_one_length(disk_radius=7, L=j)

        points = np.array(combination.calculate_circle_points())
        energy[i] = calculate_new_energy(points)
        print(energy)
        exit(0)
