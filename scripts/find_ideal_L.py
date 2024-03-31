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
    L_lst = []
    R_lst = []
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

        angles.append(round(phi_i, 9))

        # get L_i and R_i
        if i % 2 == 0:
            R_i = euclidean(next_point, current_point) / 2
            L_i = euclidean(prev_point, current_point) / 2
            L_lst.append(round(L_i,4))
            R_lst.append(round(R_i,4))
        elif i % 2 != 0:
            R_i = euclidean(prev_point, current_point) / 2
            L_i = euclidean(next_point, current_point) / 2
            L_lst.append(round(L_i,4))
            R_lst.append(round(R_i,4))
        geometry_energy += new_energy(L=L_i, phi=phi_i, R=R_i)

    print("N:", len(points)/2, "angle", set(angles))
    print("R:", set(R_lst))
    print("L:", set(L_lst))
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

    return F

def plot_energy_L(energy_lst):
    if len(energy_lst) > 100:
        plt.plot(L, energy_lst, ls='-')
    else:
        plt.plot(L, energy_lst, ls='-', marker='o', ms='3')
    plt.xlabel("L nm")
    plt.ylabel("F/N [kT]")
    # plt.grid()
    # plt.show()

if __name__ == "__main__":
    Li_default = 0.5
    parser = argparse.ArgumentParser(
        description="For each N disks finds optimal half-distance (L) for ΔN=0",
        epilog="Not all those who physics are lost"
    )
    parser.add_argument("-N", "--n_disks", required=True, type=int, help="number of proteins, REQUIRED")
    parser.add_argument("-Li", "--L_initial", type=float, default=Li_default, help="initial L for scan")
    parser.add_argument("-Lf", "--L_final", type=float, default=Li_default+5, help="final L for scan")
    parser.add_argument("-dL", "--Delta_L", type=float, default=0.5, help="dL for scan")
    parser.add_argument("-p", "--plot", action="store_true", help="plot result")
    args = parser.parse_args()

    L = np.arange(args.L_initial, args.L_final, step=args.Delta_L)
    energy_lst = np.zeros(L.shape)

    # configure geoemtry of ΔN=0
    combination = Combination(n_disks=args.n_disks)
    combination.lengths = [1]*args.n_disks # generate initial ΔN=0

    # for k in
    for i,j in enumerate(L):
        print("L", j)
        combination.modify_one_length(disk_radius=7, L=j)

        # generate geometry
        points = np.array(combination.calculate_circle_points())

        edge_color1 = 'b'
        edge_color2 = 'g'
        plot_struct = True
        if plot_struct == True:
            geom = np.vstack([points, points[0]])
            for k in range(len(geom) - 1):
                if k % 2 == 0:
                    plt.plot([geom[k][0], geom[k+1][0]], [geom[k][1], geom[k+1][1]], color=edge_color1, marker='o', ms=1, mec='r')
                else:
                    plt.plot([geom[k][0], geom[k+1][0]], [geom[k][1], geom[k+1][1]], color=edge_color2, marker='o', ms=1, mec='r')
            # plt.plot(points1[:,0], points1[:,1])
            plt.show()
            if i > 5:
                pass

        # calculate energy
        energy_lst[i] = calculate_new_energy(points)

    if args.plot:
        plot_energy_L(energy_lst)
        plt.grid()
        plt.show()
