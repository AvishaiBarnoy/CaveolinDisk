import matplotlib.pyplot as plt
import numpy as np
import os
import re

data_filename = "L_N.csv"
data = np.loadtxt(data_filename, delimiter=',')

# marker list
markers = ["o", "^", "s", "+", "h", "v", "1", "*", "8", "p", "x"]

N = [6, 8, 10, 11, 12]
Eb_k = data[:,0]
# print(Eb_k)

# for i in data:
for i,j in enumerate(data):
    plt.plot(N, data[i][1:], ls='--', marker=markers[i], label=r"$E_b/k$={KA}".format(KA=Eb_k[i]))

# plt.ylabel(r'$\tilde{F}/\epsilon$  $[\tilde{\kappa}]$', fontsize=15)
plt.xlabel('N', fontsize=15)
plt.ylabel(r'$\Delta N$', fontsize=15)
plt.title(r"$\Delta N_{min}$ with L = 2 nm with different L")

plt.legend(loc="best")
plt.xticks(N, fontsize=13)
plt.yticks(fontsize=13)
plt.grid()

# saving 
output_filename = "lowest_dn"

plt.savefig(f"{output_filename}.svg")
plt.savefig(f"{output_filename}.png")

plt.show()
