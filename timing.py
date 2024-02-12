import timeit

old = """
import numpy as np
def F_el(L, R=7, l=2):
    C = 0.4
    F = -C / (1/np.tanh(R/l) + 0.5/np.tanh(L/l))
    return F
def F_vdw(L=2, R=7, l=2):
    A = 1
    h = 2
    B = round(1 / 270 * A * h / l ** 2 * 1 / np.sqrt(R / l), 3)
    F = -B / (L / l) ** (3 / 2)
    return F
def F_tot(L, R, l=2):
    return F_el(L, R, l) + F_vdw(L, R, l)
def calc_comb_energy(length_list: list, n_disks: int, L=2, eta=0.18) -> float:
    combination_total_energy = 0.0
    for i in length_list:
        if i == 14:
            F_disk = F_tot(L=L, R=7, l=2)
            combination_total_energy += F_disk
        elif i > 14:
            F_disks = F_tot(L=L, R=7, l=2) + F_tot(L=eta, R=7, l=2) * (i - 1)
            combination_total_energy += F_disks
    return combination_total_energy

lengths = [14, 4, 14, 4, 210, 4]
calc_comb_energy(lengths, 17)
lengths = [14,14,210]
calc_comb_energy(lengths, 17)
"""
new = """
import numpy as np
def F_el(L, R=7, l=2):
    C = 0.4
    F = -C / (1/np.tanh(R/l) + 0.5/np.tanh(L/l))
    return F
def F_vdw(L=2, R=7, l=2):
    A = 1
    h = 2
    B = round(1 / 270 * A * h / l ** 2 * 1 / np.sqrt(R / l), 3)
    F = -B / (L / l) ** (3 / 2)
    return F
def F_tot(L, R, l=2):
    return F_el(L, R, l) + F_vdw(L, R, l)
def calc_comb_energy(length_list: list, n_disks: int, L=2, eta=0.18) -> float:
    combination_total_energy = 0.0
    F_disk_14 = F_tot(L=L, R=7, l=2)
    F_disk_eta = F_tot(L=eta, R=7, l=2)
    for i in length_list:
        if i == 14:
            combination_total_energy += F_disk_14
        elif i > 14:
            combination_total_energy += F_disk_14 + F_disk_eta * (i - 1)
    return combination_total_energy

lengths = [14, 4, 14, 4, 210, 4]
calc_comb_energy(lengths, 17)
lengths = [14,14,210]
calc_comb_energy(lengths, 17)
"""

n = 1000
print("new:", timeit.timeit(new, number = n))
print("old:", timeit.timeit(old, number = n))
