import matplotlib.pyplot as plt
import numpy as np
import os
import re

def calc_k(L, R=7, xi=2):
    h = 2 # nm
    a = 0.7 # nm
    k = 0.8e-19 # J
    kt = 30e-3 # N/m
    K_const = np.sqrt(k*kt/1e18) * (2/np.tanh(R/xi) + 1/np.tanh(L/xi))*1/np.tanh(L/xi) / (1/np.tanh(R/xi) + 1/np.tanh(L/xi))
    return K_const / 4.11e-21

def calc_F(L, R=7, xi=2):
    h = 2 # nm
    a = 0.7 # nm
    k = 0.8e-19 # J
    kt = 30e-3 # N/m
    depsilon = 4 # kT/nm
    A = (h/a)**2 * depsilon/np.sqrt(k*kt/1e18) * 4.11e-21
    B = 1e-3
    # print(A, B)
    F = - A / (2/np.tanh(R/xi) + 1/np.tanh(L/xi)) - B /(L/xi)**(1.5)
    return F

if __name__ == "__main__":
    L = 5
    R = 7
    xi = 2

    # data
    y = np.loadtxt("e_values.txt")
    x = np.array(range(len(y)))

    # marker list
    markers = ["o", "^", "s", "+", "h", "v", "1", "*", "8", "p", "x"]

    # different values for kappa_tilde
    L_min_lst = [0.005, 0.01, 0.015, 0.05]
    k = calc_k(L=L, R=R, xi=xi)

    F = calc_F(L=L, R=R)

    # plot
    for i,j in enumerate(L_min_lst):
        dF = calc_F(L=j) - F
        plt.plot(x, y*k + dF*x, ls='--', marker=markers[i], label=r"$\tilde{{\kappa}}={KA}$".format(KA=j))

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

        plt.savefig(f"{output_filename}.svg")
        plt.savefig(f"{output_filename}.png")

    plt.show()
