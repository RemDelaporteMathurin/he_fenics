import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import matplotx

x_right = 500  # nm
alphas = [0.05, 0.1, 0.2, 0.3, 0.4, 1, 2]

cmap = cm.get_cmap("Greens")
norm = LogNorm(vmin=1e-2, vmax=2e0)

for alpha in alphas:
    data = np.genfromtxt("k_burst_{:.2f}/r-ib.txt".format(alpha), delimiter=" ")

    x = data[:, 0] * 1e9
    ib = data[:, 1]
    plt.plot(
        x[np.where(x < x_right)][::-1],
        ib[np.where(x < x_right)][::-1],
        label="{:.2f}".format(alpha) + "$ \; k_\mathrm{burst}$",
        c=cmap(norm(alpha)),
    )

# no bursting

data = np.genfromtxt("k_burst_{:.2f}/r-ib.txt".format(0), delimiter=" ")

x = data[:, 0] * 1e9
ib = data[:, 1]
plt.plot(
    x[np.where(x < x_right)][::-1],
    ib[np.where(x < x_right)][::-1],
    label="No bursting",
    color="tab:red",
)

# plt.ylim(bottom=1e5, top=1e9)
plt.xlim(left=0, right=x_right)
plt.yscale("log")
plt.xlabel("Depth (nm)")
plt.ylabel(r"He content $\langle i_b \rangle $")
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
matplotx.line_labels()
plt.tight_layout()
plt.savefig("he_content_profile.svg")
plt.show()
