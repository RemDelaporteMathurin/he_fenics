import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
import matplotx
# try:
#     plt.rc('text', usetex=True)
#     plt.rc('font', family='serif', size=12)
# except:
#     pass


R1 = 3e-10


def radius_pure_helium(i):
    a_0 = 0.318e-9
    pi = np.pi
    r = R1
    if i > 1:
        r += (3/(4*pi)*a_0**3/10*i)**(1/3)
        r -= (3/(4*pi)*a_0**3/10)**(1/3)
    return r


def radius(i):
    a_0 = 0.318e-9
    pi = np.pi
    return ((3**0.5)/4*a_0) + \
        (3/(4*pi) * (a_0**3)/2 * abs(i)/4)**(1/3) - \
        (3/(4*pi) * (a_0**3)/2)**(1/3)


times = [10, 1, 0.1]
cmap = cm.Blues
norm = LogNorm(1e-3, 1e1)
fig, axs = plt.subplots(2, 1, sharex=True)#, figsize=(5.4, 3))

# average content
plt.sca(axs[0])
for t in times:
    data_concentration = np.genfromtxt("profiles/t={}s.csv".format(t), delimiter=",", names=True)
    data_ib = np.genfromtxt("i/t={}s.csv".format(t), delimiter=",", names=True)
    average_rad = (
        data_concentration["1"]*1 + sum([data_concentration[str(i)]*i for i in range(2, 7)]) +
        data_concentration["cb"]*data_ib["ib"]
        )[1:]
    sum_of_c_i = (data_concentration["cb"] + sum([data_concentration[str(i)] for i in range(1, 7)]))[1:]
    average_rad /= sum_of_c_i
    plt.plot(data_concentration["arc_length"][1:]*1e9, average_rad, color=cmap(norm(t)), label="{} s".format(t))
plt.yscale("log")
plt.ylabel("Average helium \n " + r"content $\langle i \rangle$ (He)")
matplotx.line_labels()
# plt.xlabel("x (m)")
# plt.legend()
# plt.show()

# average radius
plt.sca(axs[1])
for t in times:
    data_concentration = np.genfromtxt("profiles/t={}s.csv".format(t), delimiter=",", names=True)
    data_ib = np.genfromtxt("i/t={}s.csv".format(t), delimiter=",", names=True)
    average_rad = (
        data_concentration["1"]*R1 + sum([data_concentration[str(i)]*radius_pure_helium(i) for i in range(2, 7)]) +
        data_concentration["cb"]*radius(data_ib["ib"])
        )[1:]
    sum_of_c_i = (data_concentration["cb"] + sum([data_concentration[str(i)] for i in range(1, 7)]))[1:]
    average_rad /= sum_of_c_i
    plt.plot(data_concentration["arc_length"][1:]*1e9, average_rad*1e9, color=cmap(norm(t)), label="{} s".format(t))
# plt.yscale("log")
plt.ylabel("Average radius \n " + r"$\langle r \rangle$ (nm)")
plt.xlabel("x (nm)")
plt.tight_layout()
matplotx.line_labels()

for ax in axs:
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
plt.show()
