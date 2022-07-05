import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from labellines import labelLines

def radius(i):
    a_0 = 0.318e-9
    pi = np.pi
    return ((3**0.5)/4*a_0) + \
        (3/(4*pi) * (a_0**3)/2 * i/4)**(1/3) - \
        (3/(4*pi) * (a_0**3)/2)**(1/3)


folder = "i"
times = [10, 1, 0.1]

for t in times:
    data = np.genfromtxt(folder + "/t={}s.csv".format(t), delimiter=",", names=True)
    plt.plot(data["arc_length"]*1e9, radius(data["ib"])*1e9, color="black", label="{} s".format(t))

plt.xlim(5e-10*1e9, 7e-7*1e9)
x = np.logspace(np.log10(3e-10), np.log10(1e-7))
i = sp.Symbol("i")
max_i = np.array([float(sp.solve(radius(i) - val_x)[0]) for val_x in x])
plt.fill_between(x*1e9, x*1e9, np.zeros(x.shape) + x[-1]*1e9, color='grey', alpha=0.3)
plt.text(1, 49, r"$\langle r_b \rangle > x$" "\n bursting zone", fontsize=13)

plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.xlabel("x (nm)")
plt.xscale("log")
plt.ylabel(r"$\langle r_b \rangle$ (nm)")
# plt.legend()
plt.yscale("log")
plt.ylim(top=1e2)
labelLines(plt.gca().get_lines(), align=True, fontsize=14, backgroundcolor="white")
plt.show()
