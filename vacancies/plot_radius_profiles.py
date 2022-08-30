import numpy as np
import matplotlib.pyplot as plt
import matplotx


k_burst_0 = 2e3

data_bursting = np.genfromtxt(
    "k_burst_0={:.1e}/radius.txt".format(k_burst_0), delimiter=" "
)

data_no_bursting = np.genfromtxt("k_burst_0={:.1e}/radius.txt".format(0), delimiter=" ")

# remove x>xmax
x_max = 500  # nm

indexes = np.where(data_no_bursting[:, 0] < x_max * 1e-9)[0]
data_no_bursting = data_no_bursting[indexes[::-1]]

indexes = np.where(data_bursting[:, 0] < x_max * 1e-9)[0]
data_bursting = data_bursting[indexes[::-1]]

# plot

plt.plot(
    data_bursting[:, 0] * 1e9,
    data_bursting[:, 1] * 1e9,
    label="Bursting",
    color="tab:green",
)
plt.plot(
    data_no_bursting[:, 0] * 1e9,
    data_no_bursting[:, 1] * 1e9,
    label="No bursting",
    color="tab:red",
)

x = np.linspace(0.5, 1e2)
plt.fill_between(x, x, np.zeros(x.shape) + 1e2, color="grey", alpha=0.3)
plt.text(1, 12, r"$\langle r_b \rangle > x$" + "\n 'bursting zone'", fontsize=13)

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Depth (nm)")
plt.ylabel("Radius (nm)")

plt.xlim(left=1e0, right=x_max)
plt.ylim(bottom=1e0, top=1e2)
matplotx.line_labels()

plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)

plt.tight_layout()
plt.show()
