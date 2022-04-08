import re
import os
from os import path

import numpy as np
import csv

from scipy.interpolate import interp1d
from scipy.stats import linregress

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker

from inference.gp_tools import GpRegressor
from inference.gp_tools import RationalQuadratic, SquaredExponential


def fit_powerlaw(x, y):
    slope, intercept, r_value, p_value, std_err = \
        linregress(np.log10(x), np.log10(y))
    a = 10**intercept
    b = slope
    return a, b


R_p = 1.5e-9
R1 = 3e-10

def D(T):
    k_B = 8.6e-5
    return 2.95e-8*np.exp(-0.13/k_B/T)


points = []
data = []

strings = []
folder = "parametric_study"
for f in os.listdir(folder):
    strings.append(os.path.join("/" + folder, f))

count = 0
for s in strings:
    match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
    e = re.findall(match_number, s.replace(folder, ''))
    point = [float(e[i])*10**float(e[i+1]) for i in [0, 2]]
    filename = folder + "/T={:.3e}_flux={:.2e}/derived_quantities.csv".format(point[0], point[1])

    if os.path.exists(filename):
        points.append(point)
        data.append({})
        T = points[-1][0]
        data[-1]["T"] = T
        flux = points[-1][1]
        data[-1]["flux"] = flux
        point_data = np.genfromtxt(filename, delimiter=',', names=True)
        # print(point_data.dtype.names)
        c_s = flux*R_p/D(T)
        K_1_1 = 4*np.pi*(2*R1)*(2*D(T))
        l_1 = 1/(4*np.pi*R1*c_s)**0.5
        data[-1]["t"] = point_data['ts']#*c_s*K_1_1
        data[-1]["inventory"] = point_data['inventoryHem2']#/c_s/l_1
        data[-1]["inventory_interp"] = interp1d(point_data['ts'], point_data['inventoryHem2'])
        data[-1]["surface_flux"] = point_data['flux_surface_leftHem2s']#/flux
        data[-1]["surface_flux_interp"] = interp1d(point_data['ts'], point_data['flux_surface_leftHem2s'])

        data[-1]["max_ib"] = point_data['max_ibHe']
        data[-1]["max_ib_interp"] = interp1d(point_data['ts'], point_data['max_ibHe'])

        data[-1]["x_max_ib"] = point_data['x_max_ibm']#/R_p
        data[-1]["x_max_ib_interp"] = interp1d(point_data['ts'], point_data['x_max_ibm'])

        data[-1]["mean_ib"] = point_data['mean_ibHe']
        data[-1]["mean_ib_interp"] = interp1d(point_data['ts'], point_data['mean_ibHe'])

        data[-1]["total_bubbles"] = point_data['total_bubbles']#/c_s/l_1
        data[-1]["total_bubbles_interp"] = interp1d(point_data['ts'], point_data['total_bubbles'])

        if data[-1]["inventory"][-1] == data[-1]["inventory"][-4]:
            print("Warning: absolute tolerance is too low for some points")
            print(point)

if __name__ == "__main__":
    # scatter inventory
    plt.figure()
    plt.title("Helium inventory")
    xyz = np.array([[d["T"], d["flux"], d["inventory"][-1]] for d in data])
    scat = plt.scatter(xyz[:, 0], xyz[:, 1], c=xyz[:, 2])
    scat.set_norm(matplotlib.colors.LogNorm())
    plt.yscale('log')
    plt.colorbar(label="He/m2")
    plt.xlabel("T (K)")
    plt.ylabel(r"He implanted flux (He m$^{-2}$ s$^{-1}$)")

    # scatter flux
    plt.figure()
    plt.title("Outgassing flux")
    xyz = np.array([[d["T"], d["flux"], abs(d["surface_flux"][-1])] for d in data])
    scat = plt.scatter(xyz[:, 0], xyz[:, 1], c=xyz[:, 2])
    scat.set_norm(matplotlib.colors.LogNorm())
    plt.yscale('log')
    plt.colorbar(label="He/m2/s")
    plt.xlabel("T (K)")
    plt.ylabel(r"He implanted flux (He m$^{-2}$ s$^{-1}$)")

    # scatter max_ib
    plt.figure()
    plt.title("max($<i_b>$)")
    xyz = np.array([[d["T"], d["flux"], d["max_ib"][-1]] for d in data])
    scat = plt.scatter(xyz[:, 0], xyz[:, 1], c=xyz[:, 2])
    scat.set_norm(matplotlib.colors.LogNorm())
    plt.yscale('log')
    plt.colorbar(label="He/cluster")
    plt.xlabel("T (K)")
    plt.ylabel(r"He implanted flux (He m$^{-2}$ s$^{-1}$)")

    # scatter max_ib
    plt.figure()
    plt.title("max($<i_b>$) location")
    xyz = np.array([[d["T"], d["flux"], d["x_max_ib"][-1]] for d in data])
    scat = plt.scatter(xyz[:, 0], xyz[:, 1], c=xyz[:, 2]*1e9)

    plt.yscale('log')
    plt.colorbar(label="nm")
    plt.xlabel("T (K)")
    plt.ylabel(r"He implanted flux (He m$^{-2}$ s$^{-1}$)")

    plt.show()
