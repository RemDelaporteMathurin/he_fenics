import matplotlib.pyplot as plt
import csv
import numpy as np


class Cluster():
    def __init__(self, n_V, n_H, n_He):
        self.n_He = n_He
        self.n_H = n_H
        self.n_V = n_V
        self.E_H, self.E_He, self.E_Vac = None, None, None


def find_cluster(clusters, n_V, n_H, n_He):
    for cluster in clusters:
        if (cluster.n_V, cluster.n_H, cluster.n_He) == (n_V, n_H, n_He):
            return cluster
    print("Couldn't find cluster {}V.{}H.{}He".format(n_V, n_H, n_He))


def plot_binding_energies(max_He, max_H, n_Vac, binding_energy='E_H', show=True):

    x, y, z, colors = [], [], [], []
    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for cluster in clusters:
        if cluster.n_H <= max_H and cluster.n_He <= max_He and cluster.n_V == n_Vac:
            if vars(cluster)[binding_energy] is not None:
                x.append(cluster.n_He)
                y.append(cluster.n_H)
                z.append(vars(cluster)[binding_energy])
                color = cycle[cluster.n_H]
                colors.append(color)
    z = [0 if v is None else v for v in z]

    bottom = np.zeros_like(z)
    width = depth = 0.5
    fig = plt.figure()
    ax1 = fig.add_subplot(projection='3d')
    ax1.bar3d(np.array(x) - width/2, np.array(y) - depth/2, bottom, width, depth, z, color=colors, shade=True)
    plt.xlabel(r"$n_{He}$")
    plt.ylabel("$n_H$")
    if binding_energy == 'E_He':
        label = r'$E_{b_{He}}$'
    else:
        label = r'$E_{b_{H}}$'
    ax1.set_zlabel(label + ' (eV)')
    ax1.view_init(elev=20, azim=45)
    if show:
        plt.show()


clusters = []
with open('binding_energies.csv') as csvfile:
    reader = csv.reader(csvfile, delimiter=";")
    next(reader)
    for row in reader:
        n_Vac = int(row[0])
        n_He = int(row[1])
        n_H = int(row[2])
        my_cluster = Cluster(n_Vac, n_H, n_He)
        clusters.append(my_cluster)
        for i, E in enumerate(['E_H', 'E_He', 'E_Vac']):
            if row[3 + i] not in ["None", ""]:
                setattr(my_cluster, E, float(row[3 + i]))

if __name__ == "__main__":
    plot_binding_energies(6, 6, 1, 'E_He', show=False)
    plot_binding_energies(6, 6, 1, 'E_H')
