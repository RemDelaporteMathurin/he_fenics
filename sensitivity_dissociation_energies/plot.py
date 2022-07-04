import numpy as np
import matplotlib.pyplot as plt

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

def plot_concentration(i, T):
    
    lines = []
    for offset, style in zip([-0.5, 0, 0.5], ["-.", 'solid', 'dashed']):
        folder = "{}K/{}eV".format(T, offset)
        if i == "cb":
            colour = "black"
            data = np.genfromtxt(folder + "/profile_cb.csv", delimiter=",", names=True)
            label = "$c_b$"
            y = data[i]

        else:
            colour = cycle[i]
            data = np.genfromtxt(folder + "/profiles.csv", delimiter=",", names=True)
            label = "He" + str(i+1)
            y = data[str(i+1)]

        x = data["arc_length"]*1e9
        
        line, = plt.plot(x, y, linestyle=style, color=colour)
        lines.append(y)
        if offset == 0:
            line.set_label(label)
    
    plt.fill_between(x, lines[0], lines[2], color=colour, alpha=0.2)

def plot_profiles(T):
    for i in range(6):
        plot_concentration(i, T)

    # plot cb
    plot_concentration("cb", T)


fig, axs = plt.subplots(7, 3, sharex=True, sharey="row", figsize=(8, 3*3))
for axs_row, i in zip(axs, [0, 1, 2, 3, 4, 5, "cb"]):
    if i != "cb":
        label = "$c_{\mathrm{He}_"+  "{}".format(i+1) + "}$"
        if i==0:
            label += " (m$^{-3}$)"
        axs_row[0].set_ylabel(label)
    else:
        axs_row[0].set_ylabel("$c_b$")
    
    for ax, T in zip(axs_row, [500, 1000, 1500]):
        plt.sca(ax)
        plot_concentration(i, T)
        plt.yscale("log")
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
for ax in axs[-1]:
    ax.set_xlabel("x (nm)")

# fig.suptitle("Solid: +0eV ; Dashed: +0.5eV ; Dashed-point: -0.5eV")

axs[0][0].set_title("500 K")
axs[0][1].set_title("1000 K")
axs[0][2].set_title("1500 K")
plt.tight_layout()
plt.show()
