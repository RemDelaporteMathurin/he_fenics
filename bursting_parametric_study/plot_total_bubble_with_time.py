import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import matplotx

alphas = [0.05, 0.1, 0.2, 0.3, 0.4, 1, 2]
cmap = cm.get_cmap("Greens")
norm = LogNorm(vmin=1e-2, vmax=2e0)

for alpha in alphas:
    data = np.genfromtxt(
        "k_burst_{:.2f}/r-derived_quantities.csv".format(alpha),
        delimiter=",",
        names=True,
    )

    t = data["ts"]
    max_ib = data["total_bubbles"]
    plt.plot(
        t,
        max_ib,
        label="{:.2f}".format(alpha) + "$ \; k_\mathrm{burst}$",
        c=cmap(norm(alpha)),
    )
