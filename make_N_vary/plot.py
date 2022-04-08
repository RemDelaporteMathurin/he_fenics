import re
import os
from os import path

import numpy as np
import matplotlib.pyplot as plt


# find all folders
strings = []
datasets = []
folder = "."
for f in os.listdir(folder):
    strings.append(os.path.join("/" + folder, f))
for s in strings:
    match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
    e = re.findall(match_number, s)
    if len(e) != 0:
        datasets.append(int(e[0]))

# plot
fig, axs = plt.subplots(3, 1, sharex=True, figsize=[6.4, 5])
fields = ["inventoryHem2", "mean_ibHe", "total_bubbles"]
y_labels = ["He inventory  \n (He m$^{-2}$)", "mean He content in \n immobile clusters", "Bubbles inventory \n" r"(bubbles m$^{-2}$)"]
for field, y_label, ax in zip(fields, y_labels, axs):
    for nb_clusters in sorted(datasets):
        filename = "N={}/derived_quantities.csv".format(nb_clusters)
        point_data = np.genfromtxt(filename, delimiter=',', names=True)
        ax.plot(point_data['ts'], point_data[field], label="N={}".format(nb_clusters))
    ax.set_ylabel(y_label)
    ax.set_yscale("log")
    ax.set_xscale("log")
axs[-1].legend()
axs[-1].set_xlabel('Time (s)')

plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.savefig("../Figures/varying_N.svg")
plt.savefig("../Figures/varying_N.pdf")
plt.show()
