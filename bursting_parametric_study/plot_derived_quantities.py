import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import matplotx


def plot_derived_quantity(
    quantity: str,
    alphas: list,
    cmap=cm.get_cmap("Greens"),
    norm=LogNorm(vmin=1e-2, vmax=2e0),
):

    for alpha in alphas:
        data = np.genfromtxt(
            "k_burst_{:.2f}/r-derived_quantities.csv".format(alpha),
            delimiter=",",
            names=True,
        )

        t = data["ts"]
        max_ib = data[quantity]
        plt.plot(
            t,
            max_ib,
            label="{:.2f}".format(alpha) + "$ \; k_\mathrm{burst}$",
            c=cmap(norm(alpha)),
        )

    # no bursting

    data = np.genfromtxt(
        "k_burst_{:.2f}/r-derived_quantities.csv".format(0),
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


alphas = [0.05, 0.1, 0.2, 0.3, 0.4, 1, 2]

# He inventory
plot_derived_quantity("inventoryHem2", alphas=alphas)
plt.ylim(bottom=1e13)
matplotx.line_labels()
plt.ylabel(r"He inventory (m $^{-2}$ )")
plt.tight_layout()
plt.savefig("he_inventory.svg")
plt.show()

# max ib
plot_derived_quantity("max_ibHe", alphas=alphas)
# plt.ylim(bottom=1e18)
matplotx.line_labels()
plt.ylabel(r"$\max (\langle i_b \rangle )$")
plt.tight_layout()
plt.savefig("max_ibHe.svg")
plt.show()

# total bubbles
plot_derived_quantity("total_bubbles", alphas=alphas)
plt.yscale("linear")
plt.ylim(bottom=0)
matplotx.line_labels()
plt.ylabel(r"Bubbles inventory (bubbles/ m $^{-2}$ )")
plt.tight_layout()
plt.savefig("total_bubbles.svg")
plt.show()
