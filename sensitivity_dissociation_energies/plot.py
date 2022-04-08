import numpy as np
import matplotlib.pyplot as plt

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']


def plot_profiles(T):
    for i in range(6):
        lines = []
        for offset, style in zip([-0.5, 0, 0.5], ["-.", 'solid', 'dashed']):
            folder = "{}K/{}eV".format(T, offset)
            data = np.genfromtxt(folder + "/profiles.csv", delimiter=",", names=True)
            data_cb = np.genfromtxt(folder + "/profile_cb.csv", delimiter=",", names=True)
            x = data["arc_length"]
            y = data[str(i+1)]
            line, = plt.plot(x, y, linestyle=style, color=cycle[i])
            lines.append(y)
            if offset == 0:
                line.set_label("He" + str(i+1))
        plt.fill_between(x, lines[0], lines[2], color=cycle[i], alpha=0.2)

    # plot cb
    lines = []
    for offset, style in zip([-0.5, 0, 0.5], ["-.", 'solid', 'dashed']):

        folder = "{}K/{}eV".format(T, offset)
        data_cb = np.genfromtxt(folder + "/profile_cb.csv", delimiter=",", names=True)
        x = data["arc_length"]

        line, = plt.plot(data_cb["arc_length"], data_cb["cb"], color="black", linestyle=style, linewidth=1)
        lines.append(data_cb["cb"])
        if offset == 0:
            line.set_label("$C_b$")
    plt.fill_between(x, lines[0], lines[2], color="black", alpha=0.2)


fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(8, 3))
for ax, T in zip(axs, [500, 1000, 1500]):
    plt.sca(ax)
    plot_profiles(T)
    plt.xlabel("x (m)")
    plt.text(0, 1e23, "{} K".format(T))
axs[0].set_ylabel("Concentration (m$^{-3}$)")
plt.ylim(bottom=1e6, top=5e24)

fig.suptitle("Solid: +0eV ; Dashed: +0.5eV ; Dashed-point: -0.5eV")
plt.yscale("log")
axs[-1].legend(bbox_to_anchor=(1.5, 0.5), loc="center right")
# axs[0].legend(loc="lower center", ncol=4)
plt.tight_layout()
plt.show()
