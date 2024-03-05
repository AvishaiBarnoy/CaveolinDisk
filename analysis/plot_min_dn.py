import matplotlib.pyplot as plt
import numpy as np
import os
import re

data_filename = "dn_min.csv"
data = np.loadtxt(data_filename, delimiter=',')

# marker list
markers = ["o", "^", "s", "+", "h", "v", "1", "*", "8", "p", "x"]

L = [1, 2, 3, 5, 7, 10]
kappa_tilde = data[:,0]

# TODO: remove pandas
# for i in data:
for i,j in enumerate(data):
    plt.plot(L, data[i][1:], ls='--', marker=markers[i], label=r"$L_{{close}}$={KA}".format(KA=kappa_tilde[i]))

# plt.ylabel(r'$\tilde{F}/\epsilon$  $[\tilde{\kappa}]$', fontsize=15)
plt.xlabel('L [nm]', fontsize=15)
plt.ylabel(r'$\Delta N$', fontsize=15)
plt.title(r"$\Delta N_{min}$ for diffeerent $L_{close}$")

plt.legend(loc="upper right")
plt.xticks(L, fontsize=13)
plt.yticks(fontsize=13)
plt.grid()

# saving 
output_filename = "lowest_dn"

plt.savefig(f"{output_filename}.svg")
plt.savefig(f"{output_filename}.png")

plt.show()
