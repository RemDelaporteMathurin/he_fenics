import numpy as np
import matplotlib.pyplot as plt
import matplotx


def plot_derived_quantity(
    quantity: str,
):

    data = np.genfromtxt(
        "k_burst_0=2.0e+03/r-derived_quantities.csv",
        delimiter=",",
        names=True,
    )

    t = data["ts"]
    max_ib = data[quantity]
    plt.plot(
        t,
        max_ib,
        label="bursting",
    )

    # no bursting

    data = np.genfromtxt(
        "k_burst_0=0.0e+00/r-derived_quantities.csv",
        delimiter=",",
        names=True,
    )

    t = data["ts"]
    max_ib = data[quantity]
    plt.plot(
        t,
        max_ib,
        label="No bursting",
        color="tab:red",
    )

    plt.xlabel("Time (s)")

    plt.yscale("log")
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["top"].set_visible(False)


# He inventory
plot_derived_quantity("inventoryHem2")
plt.ylim(bottom=1e13)
matplotx.line_labels()
plt.ylabel(r"He inventory (m $^{-2}$ )")
plt.tight_layout()
plt.savefig("he_inventory.svg")
plt.show()

# max ib and m
fig, axs = plt.subplots(2, 1, sharex=True)
plt.sca(axs[0])
plot_derived_quantity("max_ibHe")
# plt.ylim(bottom=1e18)
matplotx.line_labels()
plt.xlabel("")
plt.ylabel(r"$\max (\langle i_b \rangle )$")

plt.sca(axs[1])
plot_derived_quantity("max_mV")
plt.ylabel(r"$\max (\langle m \rangle )$")
matplotx.line_labels()

plt.tight_layout()
plt.savefig("max_ibHe.svg")
plt.show()

# total bubbles
plot_derived_quantity("total_bubbles")
plt.yscale("linear")
plt.ylim(bottom=0)
matplotx.line_labels()
plt.ylabel(r"Bubbles inventory (bubbles/ m $^{-2}$ )")
plt.tight_layout()
plt.savefig("total_bubbles.svg")
plt.show()

# average ib m
fig, axs = plt.subplots(3, 1, sharex=True)

plt.sca(axs[0])
plot_derived_quantity("mean_ibHe")
matplotx.line_labels()
plt.ylabel(r"mean($\langle i_b \rangle$) (He)")

plt.sca(axs[1])
plot_derived_quantity("mean_mV")
plt.ylabel(r"mean($\langle m \rangle$) (V)")

matplotx.line_labels()

plt.sca(axs[2])
plot_derived_quantity("mean_ratio_HeV")
plt.yscale("linear")
plt.ylim(bottom=0)
plt.ylabel(r"mean($\langle i_b \rangle$)/mean($\langle m \rangle$) \n (He/V)")
matplotx.line_labels()

plt.xscale("log")
plt.tight_layout()
plt.savefig("mean_ib_m.svg")
plt.show()
