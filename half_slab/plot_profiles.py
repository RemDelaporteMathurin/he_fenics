import numpy as np
import matplotlib.pyplot as plt
from labellines import *

fig, (axtop, axbottom) = plt.subplots(2, 3, sharex=True, sharey="row")
folder = "profiles"
times = [0.1, 1, 10]
for ax, t in zip(axtop, times):
    data = np.genfromtxt(folder + "/t={}s.csv".format(t), delimiter=",", names=True)
    ax.plot(data["arc_length"], data["1"])
    ax.set_title("{}s".format(t))
axtop[0].set_ylabel("$C_{\mathrm{He}_1}$ (m$^{-3}$)")

for ax, t in zip(axbottom, times):
    data = np.genfromtxt(folder + "/t={}s.csv".format(t), delimiter=",", names=True)
    ax.plot(data["arc_length"], data["cb"])
    ax.set_xlabel("x (m)")
axbottom[0].set_ylabel("$C_b$ (m$^{-3}$)")
plt.show()

fig, axs = plt.subplots(1, 2, sharex=True, sharey="row")
folder = "profiles"
times = [0.0001, 0.1, 1, 10]

for t, style in zip(times, [":", "-.", "--", "solid"]):
    data = np.genfromtxt(folder + "/t={}s.csv".format(t), delimiter=",", names=True)
    axs[0].plot(data["arc_length"]*1e9, data["1"], linestyle=style, color="black", label="{} s".format(t))
axs[0].set_ylabel("Concentration (m$^{-3}$)")
axs[0].legend()
plt.sca(axs[0])
# xvals = [2e-7, 2e-7, 1e-7, 1.5e-7]
# labelLines(plt.gca().get_lines(), zorder=2.5, xvals=xvals, align=True)

for t, style in zip(times, [":", "-.", "--", "solid"]):
    data = np.genfromtxt(folder + "/t={}s.csv".format(t), delimiter=",", names=True)
    axs[1].plot(data["arc_length"]*1e9, data["cb"], linestyle=style, color="black")

axs[0].set_xlabel("x (nm)")
axs[1].set_xlabel("x (nm)")

axs[0].text(1e-7*1e9, 2.5e21, "$\mathrm{He}_1$")
axs[1].text(2e-7*1e9, 2.5e21, "Bubbles ($c_b$)")
plt.show()
