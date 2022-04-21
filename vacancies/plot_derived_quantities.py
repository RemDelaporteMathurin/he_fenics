import numpy as np
import matplotlib.pyplot as plt
import matplotx


def plot_derived_quantity(
    quantity: str,
):

    data = np.genfromtxt(
        "r-derived_quantities.csv",
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
        "no_bursting/r-derived_quantities.csv".format(0),
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

# max ib
plot_derived_quantity("max_ibHe")
# plt.ylim(bottom=1e18)
matplotx.line_labels()
plt.ylabel(r"$\max (\langle i_b \rangle )$")
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
