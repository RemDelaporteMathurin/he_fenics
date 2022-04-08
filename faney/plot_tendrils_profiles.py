import numpy as np
import matplotlib.pyplot as plt


def plot_profiles(temperature):
    main_folder = "tendrils_blondel_energies"

    folder = main_folder + "/{}K".format(temperature)
    data = np.genfromtxt(folder + "/profiles.csv", delimiter=",", names=True)
    # data_no_dissociation = np.genfromtxt(folder + "/profiles_without_dissociation.csv", delimiter=",", names=True)
    # data_cb = np.genfromtxt(folder + "/profile_cb.csv", delimiter=",", names=True)
    data_faney = np.genfromtxt(main_folder + "/faney_profiles.csv", delimiter=",", names=True)

    x = data["arc_length"]

    style_faney = "dashed"
    style_us = "solid"

    plt.plot(data_faney["Xcb_{}".format(temperature)]*1e9, data_faney["Ycb_{}".format(temperature)], color="black", linestyle=style_faney, linewidth=2)
    plt.plot(data["arc_length"]*1e9, data["cb"], color="black", linestyle=style_us, label="$C_b$", linewidth=1)

    for i in range(6):
        line = plt.plot(x*1e9, data[str(i+1)], linestyle=style_us, label="He" + str(i+1))
        plt.plot(data_faney["XHe{}_{}".format(str(i+1), temperature)]*1e9, data_faney["YHe{}_{}".format(str(i+1), temperature)], color=line[0].get_color(), linestyle=style_faney)
        # plt.plot(data_no_dissociation["arc_length"], data_no_dissociation[str(i+1)], color=line[0].get_color(), linestyle="-.")
    # plt.title("Faney (solid) / Us (dashed)")


fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(8, 2.5))
plt.subplots_adjust(wspace=0.05)
for ax, T in zip(axs, [500, 1000, 1500]):
    plt.sca(ax)
    plot_profiles(T)
    plt.xlabel("x (nm)")
    plt.text(0, 1e23, "{} K".format(T))
axs[0].set_ylabel("Concentration (m$^{-3}$)")
plt.ylim(bottom=1e6, top=5e24)


plt.yscale("log")
axs[2].legend(bbox_to_anchor=(1.5, 0.5), loc="center right")
# axs[0].legend(loc="lower center", ncol=4)
fig.tight_layout()
plt.show()
