import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
import glob
import os
import re

def dist_cm(filename):
    coords = np.loadtxt(filename)

    center_of_mass = np.mean(coords, axis=0)
    distances = np.linalg.norm(coords - center_of_mass, axis=1)

    # Calculate average distance
    avg_distance = np.mean(distances)

    return center_of_mass, avg_distance

# load data
filenames = sorted(glob.glob("dn[0-9][0-9].txt"))
x = list(range(len(filenames)))

for i in filenames:
    print(i)
    center, avg_dist = dist_cm(i)
    print("Center of Mass:", center)
    print("Average Distance from Center of Mass:", avg_dist)

m = 2
n = 4
fig, axs = plt.subplots(m, n)

for i in range(m):
    for j in range(n):
        index = i * n + j
        if index < len(filenames):
            geom = np.loadtxt(filenames[index]).tolist()
            geom.append(geom[0])
            delta_n = index + 1
            delta_n = x[index]
            axs[i, j].set_title(r"$\Delta N$ = {}".format(delta_n))
            axs[i, j].plot(*zip(*geom), marker='o', ms=2, mec='r')
            axs[i, j].title.set_size(10)

# plt.tight_layout()

for i in range(0,m):
    for j in range(0,n):
        axs[i,j].set_ylim([-35,35])
        axs[i,j].set_xlim([-35,35])
        axs[i,j].set_aspect('equal')
        axs[i,j].xaxis.set_tick_params(labelsize=8)
        axs[i,j].yaxis.set_tick_params(labelsize=8)

for ax in axs.flat:
    ax.set(xlabel='nm', ylabel='nm')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

pattern = r"L[0-9]"
cwd = os.getcwd()
output_filename = f"{re.findall(pattern, cwd)[0]}_structures"

plt.savefig(f"{output_filename}.svg")
plt.savefig(f"{output_filename}.png")
plt.show()
