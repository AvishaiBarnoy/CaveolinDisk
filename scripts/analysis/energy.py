import matplotlib.pyplot as plt
import numpy as np
import os
import re
import argparse
from scipy.spatial.distance import euclidean
import glob

def calc_k(L, R=7, k=0.4e-19, kt=20e-3, h=2, a=0.7):
    xi = np.sqrt(k/kt) * 1e9 # J / nm
    K_const = np.sqrt(k*kt/1e18) * (2/np.tanh(R/xi) + 1/np.tanh(L/xi))*1/np.tanh(L/xi) / (1/np.tanh(R/xi) + 1/np.tanh(L/xi))
    return K_const / 4.11e-21

def calc_F(L, R=7, k=0.4e-19, kt=20e-3, h=2, a=0.7, depsilon=4):
    xi = np.sqrt(k/kt) * 1e9 # J / nm
    A = (h/a)**2 * depsilon/np.sqrt(k*kt/1e18) * 4.11e-21

    B = 1 / 270 * (A * h)/ xi ** 2 * 1 / np.sqrt(R/xi)

    F = - A / (2/np.tanh(R/xi) + 1/np.tanh(L/xi)) - B /(L/xi)**(1.5)
    return F

def calculate_new_energy(file_path):
    points = np.loadtxt(file_path)
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

        geometry_energy += new_energy(L=L_i, phi=phi_i, R=R_i) # + vdw_energy(L=L_i, R=R_i)
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

def vdw_energy(L, R, k=0.4e-19, kt=20e-3, h=2, a=0.7, depsilon=4, A=1):
    """
    calculates VdW interaction between proteins.
    depracted: changed to bond energy, not distance dependent
    """
    xi = np.sqrt(k/kt) * 1e9 # J / nm

    # B = 1e-3
    B = 1 / 270 * (A * h)/ xi ** 2 * 1 / np.sqrt(R/xi)
    vdw = - B / (L/xi) ** 1.5
    return vdw

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="plot energy profiles for lowest energy geometry for each ΔN",
        epilog="Not all those who physics are lost"
    )
    parser.add_argument('-e', "--energy_method", type=str, choices=['old', 'new'], default="new", help="""old assumes that in
initial geometry the angles areclose to ideal angle and in the new one they aren't, the methods use different energy functions
albeit a little similar.""")
    parser.add_argument('-s', '--silent', action="store_true", help="supresses plt.show()")
    args = parser.parse_args()

    L = 5
    R = 7

    # data
    y = np.loadtxt("e_values.txt")
    # TODO: use regex to extract dN 
    x = np.array(range(len(y)))

    # marker list
    markers = ["o", "^", "s", "+", "h", "v", "1", "*", "8", "p", "x"]

    # different values for kappa_tilde
    L_min_lst = [0.003, 0.005, 0.01, 0.015, 0.05]

    if args.energy_method == 'old':
        k = calc_k(L=L, R=R)
        F = calc_F(L=L, R=R)

        # plot
        for i,j in enumerate(L_min_lst):
            dF = calc_F(L=j) - F
            plt.plot(x, y*k + dF*x, ls='--', marker=markers[i], label=r"${{\kappa}}={KA}$".format(KA=j))

    elif args.energy_method == 'new':
        L_min_lst = [0.003, 0.004, 0.0045, 0.005]
        F_bind = [8, 9, 10, 11, 12, 13, 14]
        F_bind = range(5, 16, 2)
        # load data
        filenames = sorted(glob.glob("dn[0-9][0-9].txt"))
        x = list(range(len(filenames)))
        print(filenames)
        dn0_path = "dn00.txt"
        F_ref = calculate_new_energy(dn0_path)
        print(f"reference state (ΔN=0) energy: {round(F_ref,3)} kT")
        x = np.array(x)
        # for i,j in enumerate(L_min_lst):
        for i,j in enumerate(F_bind):
            dF = np.array([calculate_new_energy(filename) for filename in filenames]) - F_ref
            # print("dF", dF)
            # F_tot = dF + x * vdw_energy(L=j, R=7)
            # print("L_min:", j)
            # print("de:", x * vdw_energy(L=j, R=7), '\n')

            k = 0.4e-19
            kt = 20e-3
            K = np.sqrt(k*kt) * 1e-9 / 4.11e-21 # kT/nm
            F_tot = dF - x * F_bind[i]

            plt.plot(x, F_tot, ls='--', marker=markers[i], label=f"{j}")

    plt.ylabel(r'$F$', fontsize=15)
    plt.xlabel(r'$\Delta N$', fontsize=15)

    plt.legend()
    plt.xticks(x, fontsize=13)
    plt.yticks(fontsize=13)
    plt.grid()

    save = True
    if save == True:
        # unique output name
        pattern = r"L[0-9]"
        cwd = os.getcwd()
        output_filename = f"{re.findall(pattern, cwd)[0]}_energy"

        # plt.savefig(f"{output_filename}.svg")
        plt.savefig(f"{output_filename}.png")

    if not args.silent:
        plt.show()
