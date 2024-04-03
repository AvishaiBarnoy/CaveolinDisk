import numpy as np
import matplotlib.pyplot as plt

def calc_ideal_angle(L, R=7, h=2, a=0.7, depsilon=4, k=0.4e-19, kt=20e-3):
    xi = np.sqrt(k/kt) * 1e9 # J / nm
    f_param = h/a * depsilon * 1/np.sqrt(k*kt/1e18) * 4.11e-21
    ideal_angle = np.pi- f_param * 1 / (2/np.tanh(R/xi) + 1/np.tanh(L/xi))
    return ideal_angle


if __name__ == "__main__":
    L = np.linspace(0.5, 2, 100)
    R = [7, 14, 21, 28, 35]

    for i in R:
        angles = calc_ideal_angle(L=L, R=i)
        plt.plot(L, angles, lw=2, label=f"R={i}nm")
    # plt.plot(L, calc_ideal_angle(L=L, R=7), lw=2, label=f"R=7nm")
    # plt.plot(L, calc_ideal_angle(L=L, R=14), lw=2, label=f"R=14nm")
    # plt.plot(L, calc_ideal_angle(L=L, R=21), lw=2, label=f"R=21nm")

    plt.legend(fontsize=11)
    plt.xlabel(r"$\phi^*$", fontsize=13)
    plt.ylabel("L [nm])", fontsize=13)
    plt.grid()
    plt.show()
