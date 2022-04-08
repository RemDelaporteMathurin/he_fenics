import re
import os
from os import path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import csv
from scipy.interpolate import griddata

radius = 1e-9
points = []
data = []

# extract high temp data
strings = []
folder = "bubble_size"
for f in os.listdir(folder):
    strings.append(os.path.join(folder, f))
print("Number of simulation points: {:.0f}".format(len(strings)))
count = 0
for s in strings:
    match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
    e = re.findall(match_number, s)
    time = None
    with open(s, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            if float(row[1]) >= radius and float(row[0]) > 1e-4:
                time = float(row[0])
                break
    if time is not None:
        points.append([float(e[i])*10**float(e[i+1]) for i in [0, 2]]+[time])
        data.append({})
        data[-1]["T"] = points[-1][0]
        data[-1]["flux"] = points[-1][1]
        data[-1]["time"] = points[-1][2]

points = np.array(points)

Ts = sorted(np.unique(points[:, 0]))
fluxes = sorted(np.unique(points[:, 1]))

for T in Ts:
    y = []
    x = []
    for f in fluxes:
        count = 0
        for d in data:
            count += 1
            if d["T"] == T and d["flux"] == f:
                y.append(d["time"])
                x.append(d["flux"])
                break
    if len(y) > 2:
        plt.plot(x, y, marker='o', label="T={:.0f} K".format(T))
plt.title("Time to get {:.0e} m radius bubbles".format(radius))
plt.xlabel("flux (He/m2/s)")
plt.ylabel("time (s)")
plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.minorticks_on()
plt.grid(which='minor', alpha=0.3)
plt.grid(which='major', alpha=0.7)
plt.show()

Ts = sorted(np.unique(points[:, 0]))
fluxes = sorted(np.unique(points[:, 1]))
for f in fluxes:
    y = []
    x = []
    for T in Ts:
        for d in data:
            if d["T"] == T and d["flux"] == f:
                y.append(d["time"])
                x.append(d["T"])
    if len(y) > 2:
        plt.plot(x, y, marker='o', label="flux={:.0e} He/m2/s".format(f))
plt.title("Time to get {:.0e} m radius bubbles".format(radius))
plt.xlabel("T (K)")
plt.ylabel("time (s)")
plt.legend()
plt.yscale("log")
plt.minorticks_on()
plt.grid(which='minor', alpha=0.3)
plt.grid(which='major', alpha=0.7)
plt.show()


def get_time(points):
    ''' returns a 1D array'''
    values = []
    for p in points:
        for d in data:
            if d["T"] == p[0] and d["flux"] == p[1]:
                values.append(p[2])
                break
    return np.asarray(values)

# 3D scatter
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(points[:, 0], np.log10(points[:, 1]), np.log10(points[:, 2]))
plt.show()

# create 2D contour

def scientificNotation(value):
    if value == 0:
        return '0'
    else:
        e = np.log10(np.abs(value))
        m = np.sign(value) * 10 ** (e - int(e))
        return r'${:.0f} \times 10^{{{:d}}}$'.format(m, int(e))

x = np.linspace(150, 1000, num=10)
y = np.linspace(16, 21, num=15)
X, Y = np.meshgrid(x, y)

# grid the data
grid_z = griddata(
    (points[:, 0], np.log10(points[:, 1])),
    np.log10(get_time(points)),
    (X, Y), method='cubic', fill_value=np.nan)
locator = ticker.LogLocator(base=10)

levels = np.logspace(
    np.log10(np.min(points[:, 2])),
    # np.log10(1e1),
    np.log10(np.max(points[:, 2])),
    1000)

plt.title("Time to get {:.0e} m radius bubbles".format(radius))
CF = plt.contourf(X, 10**Y, 10**grid_z, levels=levels, locator=locator)
CS = plt.contour(X, 10**Y, 10**grid_z, levels=10, locator=locator, colors="white")
CLS = plt.clabel(CS, inline=True, fontsize=10, fmt=scientificNotation)
# plt.scatter(X2, Y2, color="red", marker="+")
plt.scatter(points[:, 0], points[:, 1], marker="+", color="grey", label="Simulation points")
plt.colorbar(CF, label=r"Time (s)", ticks=locator, extend="max")
plt.yscale("log")
plt.legend()
plt.xlabel("T (K)")
plt.ylabel("flux (He/m2/s)")
plt.show()
