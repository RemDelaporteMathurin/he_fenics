import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


def radius(i):
    a_0 = 0.318e-9
    pi = np.pi
    return ((3**0.5)/4*a_0) + \
        (3/(4*pi) * (a_0**3)/2 * i/4)**(1/3) - \
        (3/(4*pi) * (a_0**3)/2)**(1/3)


data = np.genfromtxt("profile_i.csv", delimiter=",", names=True)
arc_length = data["arc_length"]*1e9
radius =  radius(data["ib"])*1e9
if __name__=="__main__":
    plt.plot(arc_length, radius, color="black")

    plt.xlim(5e-10*1e9, 7e-7*1e9)

    plt.xlabel("x (nm)")
    plt.xscale("log")
    plt.ylabel(r"$\langle r_b \rangle$ (nm)")
    plt.yscale("log")
    plt.ylim(top=1e2)
    plt.show()
