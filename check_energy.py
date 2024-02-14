import numpy as np
import matplotlib.pyplot as plt

def F_tot(L, R=7, l=2):
    A = 1  
    h = 2  
    B = 1 / 270 * A * h / l ** 2 * 1 / np.sqrt(R / l)
    B = 1e-3
    F_vdw = - B / (L*2 / l) ** (3/2)
            
    C = 0.4
    F_el = - C / (1/np.tanh(R/l) + 0.5 / np.tanh(L/l))
    return F_vdw + F_el

f_ela = lambda r: - 0.4 / (1/np.tanh(1) + 0.5 / np.tanh(r/2))
f_vdw = lambda r, B: - B / (r*2/ 2) ** (3/2)

def calc_comb_energy(lengths) -> float:
    energy = 0.0
    eta = 0.18 / 2
    L = 2

    length_list = np.array(lengths) // 14
    non_zero = np.count_nonzero(length_list)
    energy += non_zero * F_tot(L=2)

    n_close_inter = length_list[length_list != 0] - 1
    energy += n_close_inter.sum() * F_tot(L=eta)
    print(f"F_tot(eta): {F_tot(L=eta)}, eta={eta}")
    return energy

if __name__ == "__main__":
    ene_basic = F_tot(L=2)
    print(f"basic: {ene_basic}")
    
    eta = 0.09 / 10
    ene = F_tot(L=eta)
    print(f"close: {ene}, eta: {eta}") 

    # lengths = [14, 4] * 17
    # comb = calc_comb_energy(lengths)
    # print(f"comb no connect: {comb}")
    
    lengths = [14, 4] * 13 + [28, 4, 42, 4]
    lengths = [14, 4] * 15 + [28, 4]
    comb = calc_comb_energy(lengths)
    # print(lengths)
    print(f"comb with connect: {comb}")
    
    r = np.linspace(0.1, 2.1, 100)
    e_ela = f_ela(r)
    B = 1e-2
    e_vdw = f_vdw(r, B)
    plt.plot(r, e_ela, ls='--', lw=2.5, label='elastic')
    plt.plot(r, e_vdw, ls='--', lw=2.5, label=f'vdw B={B}')
    plt.plot(r, f_vdw(r, 5e-3), ls='--', label='vdw B=5e-3')
    plt.plot(r, e_vdw+e_ela, ls='-', lw=3, label='total')
    plt.grid(); plt.legend(fontsize=13)
    plt.xlabel('half distance between disks [nm]', fontsize=13)
    plt.ylabel('energy [kT]', fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.title("interaction energy between adjacent disks $F_{tot}$\n 0.1nm to 2.1nm", fontsize=18)
    plt.show()
