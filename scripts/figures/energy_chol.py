import numpy as np
import matplotlib.pyplot as plt

def energy(k=0.4e-19, kt=20e-3, zeta=-1/2.9, a=0.7, L=2, c0=0.3):
    # k = 0.4e-19 # J, bending rigidity
    # TODO: check units
    # kt = 20e-3 # J/m2, tilt modulus 
    # zeta = 1 # nm-1, spontaneous curvature
    # a = 0.7  # nm, headgroup size
    # L = 2 # nm, half-distance
    # c0 = 0.3 # mole fraction of cholesterol in reservoir

    xi = np.sqrt(k/kt) * 1e9
    rho = R/xi
    l = L/xi
    K = np.sqrt(k*kt) * 1e-9 / 4.11e-21 # kT/nm
    Ein = K * (1/np.tanh(rho) - 0.5*zeta**2 * a*c0*(1-c0) * (rho + np.sinh(rho)*np.cosh(rho))/(np.sinh(rho)**2))
    Eout = K * (1/np.tanh(l) - 0.5*zeta**2 * a*c0*(1-c0) * (l + np.sinh(l)*np.cosh(l))/(np.sinh(l)**2))

    F = Ein + Eout
    return

def plot_line():
    plt.show()

if __name__ == "__main__":
    pass
