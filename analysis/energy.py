import matplotlib.pyplot as plt
import numpy as np
import os
import re

# data
y = np.loadtxt("e_values.txt")
y = y - y[0]
x = list(range(len(y)))

# marker list
markers = ["o", "^", "s", "+", "h", "v", "1", "*", "8", "p", "x"]

# different values for kappa_tilde
scaling = [1, 2, 5, 10, 20, 30, 50]

# plots
for i,j in enumerate(scaling):
    plt.plot(y*j - x, ls='--', marker=markers[i], label=r"$\tilde{{\kappa}}={KA}$".format(KA=j))

plt.ylabel(r'$\tilde{F}/\epsilon$  $[\tilde{\kappa}]$', fontsize=15)
plt.xlabel(r'$\Delta N$', fontsize=15)

plt.legend()
plt.xticks(x, fontsize=13)
plt.yticks(fontsize=13)
plt.grid()

# unique output name
pattern = r"L[0-9]"
cwd = os.getcwd()
output_filename = f"{re.findall(pattern, cwd)[0]}_energy"

plt.savefig(f"{output_filename}.svg")
plt.savefig(f"{output_filename}.png")
plt.show()
