import numpy as np
import matplotlib.pyplot as plt
import plot_i

try:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)
except:
    pass

fig, axs = plt.subplots(2, 1, sharex=True)

# density vs depth
plt.sca(axs[0])
data = np.genfromtxt("formatted_data/output_Density2D(depth)12_files.txt", delimiter=",", names=True)
depth = data["Depth_nm"]
density = data["Density2D_1mum2"]
error_bars = data["errDensity2D_1mum2"]
m3_to_um3 = 1e18
um3_to_m3 = 1/m3_to_um3

plt.scatter(depth, density/um3_to_m3, alpha=0.5, label="Experiment")
plt.errorbar(depth, density/um3_to_m3, yerr=error_bars/um3_to_m3, linestyle="None", capsize=3, capthick=1, alpha=0.5)

data = np.genfromtxt("profile_cb.csv", delimiter=",", names=True)
plt.plot(data["arc_length"]*1e9, data["cb"], label="Model")
# plt.xlabel("Depth (nm)")
plt.ylabel("Density (m$^{-3}$)")
plt.yscale("log")
# plt.legend()

# diamater vs depth
plt.sca(axs[1])
data = np.genfromtxt("formatted_data/output_Diameter(depth)12_files.txt", delimiter=",", names=True)
depth = data["Depth_nm"]
diameter = data["Diameter_nm"]
radius = diameter/2
error_bars = data["errDiameter_nm"]

plt.scatter(depth, radius, alpha=0.5, label="Experiment")
plt.errorbar(depth, radius, yerr=error_bars/2, linestyle="None", capsize=3, capthick=1, alpha=0.5)
plt.xlabel("Depth (nm)")
plt.ylabel("Radius (nm)")
# plt.xscale("log")
plt.xlim(left=0, right=100)
plt.yscale("log")

plt.plot(plot_i.arc_length, plot_i.radius, label="Model")
plt.plot(plot_i.arc_length, 0.1*plot_i.radius, label=r"Model $\times 0.1$", linestyle="--", color="tab:orange")
plt.legend(loc="upper right")
plt.show()
