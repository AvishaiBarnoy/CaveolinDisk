import numpy as np
import matplotlib.pyplot as plt
# import scipy as sp
from scipy import optimize

def ideal_phi(k=0.4e-19, kt=20e-3, zeta=-1/3, a=0.7, L=2, c0=0.3, depsilon=4, R=7, h=2):
    xi = np.sqrt(k/kt) * 1e9
    rho = R/xi
    l = L/xi
    K = np.sqrt(k*kt) * 1e-9 / 4.11e-21 # kT/nm
    kT = 1

    Ein =   K * (1/np.tanh(rho) - 0.5 * k/4.11e-21 * zeta**2/kT * a*c0*(1-c0) * (rho + np.sinh(rho)*np.cosh(rho))/(np.sinh(rho)**2))
    Eout =  K * (1/np.tanh(l)   - 0.5 * k/4.11e-21 * zeta**2/kT * a*c0*(1-c0) * (l + np.sinh(l)*np.cosh(l))/(np.sinh(l)**2))

    phi_star = h/a * depsilon/(Ein + 0.5*Eout)
    return phi_star

def F_tot(L, k=0.4e-19, zeta=-1/2.9, phi=ideal_phi(), depsilon=4, c0=0.3, kt=20e-3, a=0.7, R=7, h=2):
    xi = np.sqrt(k/kt) * 1e9
    rho = R/xi
    l = L/xi
    K = np.sqrt(k*kt) * 1e-9 / 4.11e-21 # kT/nm
    kT = 1 # kT

    Ein =   K * (1/np.tanh(rho) - 0.5 * k/4.11e-21 * zeta**2/kT * a*c0*(1-c0) * (rho + np.sinh(rho)*np.cosh(rho))/(np.sinh(rho)**2))
    Eout =  K * (1/np.tanh(l)   - 0.5 * k/4.11e-21 * zeta**2/kT * a*c0*(1-c0) * (l + np.sinh(l)*np.cosh(l))/(np.sinh(l)**2))

    # F = 0.5 * (Ein + 0.5*Eout) * phi**2 - h/a * depsilon * phi
    F = (2*Ein+Eout)/(2*Ein+2*Eout)*Eout*phi**2 - (Eout)/(Ein+Eout)*h/a*depsilon*phi - (h/a*depsilon)**2/(2*Ein+2*Eout)
    return F

if __name__ == "__main__":
    L = np.linspace(1, 8.5, 200)
    # k_lst = [0.4e-19, 1.5e-19]
    k_lst = np.linspace(0.4e-19, 1.5e-19, num=5)
    fig, ax = plt.subplots()

    # TODO: replace with nested loops
    L_min1 = []
    k_min1 = []
    L_min2 = []
    k_min2 = []
    L_min3 = []
    k_min3 = []

    for i in k_lst:
        # print("k", i)
        phi_star = ideal_phi(k=i)
        energy = F_tot(phi=phi_star, k=i, L=L)
        ax.plot(L, energy, label=f"$\kappa$={round(i,20)} J")
        cg = optimize.minimize(F_tot, x0=0.5, args=(i), method='cg')
        ax.scatter(cg.x, F_tot(L=cg.x, k=i, phi=phi_star), s=20)

        # F_tot(L, k=0.4e-19, zeta=-1/2.9, phi=ideal_phi(), depsilon=4, c0=0.3, kt=20e-3, a=0.7, R=7, h=2)
        cg1 = optimize.minimize(F_tot, x0=0.5, args=(i,-1/2.5))
        cg2 = optimize.minimize(F_tot, x0=0.5, args=(i,-1/2))
        L_min1.append(cg.x)
        k_min1.append(i)
        L_min2.append(cg1.x)
        k_min2.append(i)
        L_min3.append(cg2.x)
        k_min3.append(i)

    ax.set_xlabel(r"$L [nm]$")
    ax.set_ylabel(r"$F_{tot}$ [$k_B T$]")
    ax.legend()
    ax.grid()
    output_name = "F_L.png"
    fig.savefig(output_name)
    plt.show()

    plt.plot(k_min1, L_min1, "o-", ms=3)
    plt.plot(k_min2, L_min2, "o-", ms=3)
    plt.plot(k_min3, L_min3, "o-", ms=3)
    zeta = [1/3, 1/2.5, 1/2]
    plt.legend([f"$\zeta$={round(z,2)} $nm^{{-1}}$]" for z in zeta])
    plt.grid()
    plt.xlabel(r"$\kappa\ [J]$")
    plt.ylabel(r"$L^*\ [nm]$")
    plt.savefig("k_L.png")
    plt.show()
    # for i in list(zip(k_min3,L_min3)):
        # print(i)

    phi = np.linspace(0, np.pi/2, num=100)
    # F_tot(L, k=0.4e-19, zeta=-1/2.9, phi=ideal_phi(), depsilon=4, c0=0.3, kt=20e-3, a=0.7, R=7, h=2)
    k1 = 1.5e-19
    plt.plot(phi, F_tot(L=1, phi=phi, k=k1), label="L=1 nm")
    plt.plot(phi, F_tot(L=2, phi=phi, k=k1), label="L=2 nm")
    plt.plot(phi, F_tot(L=5, phi=phi, k=k1), label="L=5 nm")
    plt.plot(phi, F_tot(L=7, phi=phi, k=k1), label="L=7 nm")
    plt.xlabel(r"$\varphi\ [rad]$")
    plt.ylabel(r"$F_{tot}\ [k_B T]$")
    plt.grid(); plt.legend()
    plt.savefig(r"F_vs_varphi.png")
    plt.show()
