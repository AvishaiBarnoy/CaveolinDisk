from combinations import Combination
import numpy as np
import argparse
from scipy.spatial.distance import euclidean
import matplotlib.pyplot as plt
"""
1. For each L, build geometry
2. Calculate energy for each geometry
3. move to next L
"""

def calculate_new_energy(points):
    geometry_energy = 0
    angles = []
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

        angles.append(round(phi_i, 4))

        # get L_i and R_i
        if i % 2 == 0:
            R_i = euclidean(next_point, current_point) / 2
            L_i = euclidean(prev_point, current_point) / 2
        elif i % 2 != 0:
            R_i = euclidean(prev_point, current_point) / 2
            L_i = euclidean(next_point, current_point) / 2

        geometry_energy += new_energy(L=L_i, phi=phi_i, R=R_i)

    # print("N:", len(points)/2, set(angles))
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
    # print(round(phi,4))
    f1 = 0.5 * K * l * (2*r + l) / (r + l) * phi**2
    f2 = l / (r + l) * e * phi
    f3 = 0.5 / (r + l) * e**2 / K

    F = f1 - f2 - f3

    return F

def plot_energy_L(energy_lst, N, L):
    if len(energy_lst) >= 100:
        plt.plot(N, energy_lst, ls='-')
    else:
        plt.plot(N, energy_lst, ls='-', marker='o', ms='3', label=f"L: {round(L,2)}")
    plt.xlabel("N disks", fontsize=13)
    plt.ylabel(f"F/N [kT]", fontsize=13)
    plt.tick_params(axis='both', which='minor', labelsize=11)
    # plt.title(f"F(N) for L = {L} nm")
    # plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="For each N disks finds optimal half-distance (L) for ΔN=0",
        epilog="Not all those who physics are lost"
    )
    parser.add_argument("-L", required=True, type=float, help="half-distance with dN=0, REQUIRED")
    parser.add_argument("-Ni", "--N_initial", type=int, default=4, help="initial N for scan")
    parser.add_argument("-Nf", "--N_final", type=int, default=15, help="final N for scan")
    parser.add_argument("-p", "--plot", action="store_true", help="plot result")
    args = parser.parse_args()

    N = np.arange(args.N_initial, args.N_final+1, dtype=int, step=1)
    energy_lst = np.zeros(N.shape)

    for L in np.linspace(0.5, 10, 10):
        for i,j in enumerate(N):
            # configure geoemtry of ΔN=0
            combination = Combination(n_disks=j)
            combination.lengths = [1]*j # generate initial ΔN=0
            combination.modify_one_length(disk_radius=7, L=L)

            # generate geometry
            points = np.array(combination.calculate_circle_points())

            # calculate energy
            energy_lst[i] = calculate_new_energy(points) / j # normalized

            # plot individual geometries, using zebra
            plot_struct = False
            if plot_struct == True:
                edge_color1 = 'b'
                edge_color2 = 'g'
                geom = np.vstack([points, points[0]])
                for k in range(len(geom) - 1):
                    if k % 2 == 0:
                        plt.plot([geom[k][0], geom[k+1][0]], [geom[k][1], geom[k+1][1]], color=edge_color1, marker='o', ms=1, mec='r')
                    else:
                        plt.plot([geom[k][0], geom[k+1][0]], [geom[k][1], geom[k+1][1]], color=edge_color2, marker='o', ms=1, mec='r')
                plt.show()

        plot_energy_L(energy_lst=energy_lst, N=N, L=L)

    if args.plot:
        plt.legend()
        plt.grid()
        plt.show()
