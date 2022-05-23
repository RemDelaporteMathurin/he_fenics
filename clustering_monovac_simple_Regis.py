from fenics import *
import sympy as sp
from src import meshing, initialising, solving
import numpy as np
import matplotlib.pyplot as plt

x = sp.Symbol("x[0]")
t = sp.Symbol("t")
k_B = 8.617e-5  # eV/K


# 3e16; 7e16; 2.9e17; 5e17 (He.m^-2.s^-1)
fluxflux=5e17

N = 4
nb_clusters = 2 * N + 1
# Files
files = []
folder = "kornelsen/results"
for i in range(0, nb_clusters):
    if i <= N - 1:
        m = 0
        n = i + 1
    else:
        m = 1
        n = i - N
    files.append(
        XDMFFile("clustering/" + folder + "/c_" + str(n) + "He" + str(m) + "V.xdmf")
    )
    files[i].parameters["flush_output"] = True
    files[i].parameters["rewrite_function_mesh"] = False

# Defining mesh
size = 2e-3

mesh_parameters = {
    "size": size,
    "initial_number_of_cells": 1200,
    "refinements": [
        {"x": 2e-6, "cells": 500},
        {"x": 10e-9, "cells": 200},
    ],
}
mesh = meshing.mesh_and_refine(mesh_parameters)
vm, sm = meshing.subdomains(
    mesh,
    {
        "mesh_parameters": mesh_parameters,
        "materials": [{"borders": [0, size], "id": 1}],
    },
)
dt = Constant(2)
dx = Measure("dx", domain=mesh)
# Defining functions
V = VectorFunctionSpace(mesh, "P", 1, nb_clusters)
V_DG1 = FunctionSpace(mesh, "DG", 1)
V_DG0 = FunctionSpace(mesh, "DG", 0)
c = Function(V)
v = TestFunction(V)
sols = list(split(c))
sols.insert(0, 0)
sols = np.reshape(sols, (-1, int(len(sols) / 2)))
test_func = list(split(v))
test_func.insert(0, 0)
test_func = np.reshape(test_func, (-1, int(len(test_func) / 2)))
ini = []

# Initial conditions
width = 1.1e-9
center = 2.0e-9
distribution = (
    1 / (width * (2 * 3.14) ** 0.5) * sp.exp(-0.5 * ((x - center) / width) ** 2)
)
ini = ({"value": fluxflux * distribution, "component": N},)  # laflux

c_n = initialising.initialise_solutions({"initial_conditions": ini}, V)
c_n_back = Function(V)
c_n_back.assign(c_n)
prev_sols = list(split(c_n))
prev_sols.insert(0, 0)
prev_sols = np.reshape(prev_sols, (-1, int(len(prev_sols) / 2)))
bcs = []

bcs.append(DirichletBC(V.sub(0), Constant(0), sm, 1))
bcs.append(DirichletBC(V.sub(0), Constant(0), sm, 2))
bcs.append(DirichletBC(V.sub(1), Constant(0), sm, 1))
bcs.append(DirichletBC(V.sub(1), Constant(0), sm, 2))
bcs.append(DirichletBC(V.sub(2), Constant(0), sm, 1))
bcs.append(DirichletBC(V.sub(2), Constant(0), sm, 2))
bcs.append(DirichletBC(V.sub(3), Constant(0), sm, 1))
bcs.append(DirichletBC(V.sub(3), Constant(0), sm, 2))
bcs.append(DirichletBC(V.sub(N), Constant(0), sm, 1))
bcs.append(DirichletBC(V.sub(N), Constant(0), sm, 2))
# Diffusion coeff
T = Constant(300)
diff = np.zeros((2, N + 1)).tolist()
diff[0][1] = 2.95e-8 * exp(-0.13 / k_B / T)
diff[0][2] = 3.24e-8 * exp(-0.2 / k_B / T)
diff[0][3] = 2.26e-8 * exp(-0.25 / k_B / T)
diff[0][4] = 1.68e-8 * exp(-0.2 / k_B / T)
diff[1][0] = 177e-8 * exp(-1.29 / k_B / T)  # Diffusivity of V1

# Capture radii
a_0 = 0.318e-9
r_He1 = 3e-10
r_V1 = 3**0.5 / 4 * a_0
r = [[0, r_He1], [r_V1]]

i = 1
while len(r[0]) != N + 1:
    i += 1
    r[0].append(
        r_He1
        + (3 / 4 / pi * a_0**3 / 10 * i) ** (1 / 3)
        - (3 / 4 / pi * a_0**3 / 10) ** (1 / 3)
    )
i = 1
while len(r[1]) != N + 1:
    i += 1
    r[1].append(
        r_V1
        + (3 / 4 / pi * a_0**3 / 2 * 1) ** (1 / 3)
        - (3 / 4 / pi * a_0**3 / 2) ** (1 / 3)
    )


# Reaction rates

# Dissociation energies (emission of 1 He)
E_b_He = [[None, None, 1, 1.1, 1.2, 1.3], [None, 4, 3, 2.75, 2.45, 2.3]]

density = Constant(6.3e28)
reactions = [
    # pure He cluster
    {
        "reactives": [(0, 1), (0, 1)],
        "k+": 4 * pi * (diff[0][1] + diff[0][1]) * (r[0][1] + r[0][1]),
        "k-": 4
        * pi
        * (diff[0][1] + diff[0][1])
        * (r[0][1] + r[0][1])
        * density
        * exp(-E_b_He[0][2] / k_B / T),
    },
    {
        "reactives": [(0, 1), (0, 2)],
        "k+": 4 * pi * (diff[0][1] + diff[0][2]) * (r[0][1] + r[0][2]),
        "k-": 4
        * pi
        * (diff[0][1] + diff[0][2])
        * (r[0][1] + r[0][2])
        * density
        * exp(-E_b_He[0][3] / k_B / T),
    },
    {
        "reactives": [(0, 1), (0, 3)],
        "k+": 4 * pi * (diff[0][1] + diff[0][3]) * (r[0][1] + r[0][3]),
        "k-": 4
        * pi
        * (diff[0][1] + diff[0][3])
        * (r[0][1] + r[0][3])
        * density
        * exp(-E_b_He[0][4] / k_B / T),
    },
    # absorbtion of 1 He1
    {
        "reactives": [(0, 1), (1, 0)],
        "k+": 4 * pi * (diff[0][1] + diff[1][0]) * (r[0][1] + r[1][0]),
        "k-": 4
        * pi
        * (diff[0][1] + diff[1][0])
        * (r[0][1] + r[1][0])
        * density
        * exp(-E_b_He[1][1] / k_B / T),
    },
    {
        "reactives": [(0, 1), (1, 1)],
        "k+": 4 * pi * (diff[0][1] + diff[1][1]) * (r[0][1] + r[1][1]),
        "k-": 4
        * pi
        * (diff[0][1] + diff[1][1])
        * (r[0][1] + r[1][1])
        * density
        * exp(-E_b_He[1][2] / k_B / T),
    },
    {
        "reactives": [(0, 1), (1, 2)],
        "k+": 4 * pi * (diff[0][1] + diff[1][2]) * (r[0][1] + r[1][2]),
        "k-": 4
        * pi
        * (diff[0][1] + diff[1][2])
        * (r[0][1] + r[1][2])
        * density
        * exp(-E_b_He[1][3] / k_B / T),
    },
    {
        "reactives": [(0, 1), (1, 3)],
        "k+": 4 * pi * (diff[0][1] + diff[1][3]) * (r[0][1] + r[1][3]),
        "k-": 4
        * pi
        * (diff[0][1] + diff[1][3])
        * (r[0][1] + r[1][3])
        * density
        * exp(-E_b_He[1][4] / k_B / T),
    },
    # absorption of 1 V1
    {
        "reactives": [(0, 2), (1, 0)],
        "k+": 4 * pi * (diff[0][2] + diff[1][0]) * (r[0][1] + r[1][0]),
        "k-": 0,  # 4*pi*(diff[0][2] + diff[1][0])*(r[0][1] + r[1][0])*density*exp(-E_b[1][2]/k_B/T),
    },
    {
        "reactives": [(0, 3), (1, 0)],
        "k+": 4 * pi * (diff[0][3] + diff[1][0]) * (r[0][1] + r[1][0]),
        "k-": 0,  # 4*pi*(diff[0][3] + diff[1][0])*(r[0][1] + r[1][0])*density*exp(-E_b[1][3]/k_B/T),
    },
    {
        "reactives": [(0, 4), (1, 0)],
        "k+": 4 * pi * (diff[0][4] + diff[1][0]) * (r[0][1] + r[1][0]),
        "k-": 0,  # 4*pi*(diff[0][4] + diff[1][0])*(r[0][1] + r[1][0])*density*exp(-E_b[1][4]/k_B/T),
    },
]
F = 0
center = 2.0e-9
width = 1.1e-9
distribution = (
    1 / (width * (2 * 3.14) ** 0.5) * sp.exp(-0.5 * ((x - center) / width) ** 2)
)
flux_He = fluxflux # laflux
source = Expression(
    "t < implantation_time ?" + sp.printing.ccode(distribution) + "*flux_He: 0",
    degree=2,
    t=0,
    flux_He=flux_He,
    implantation_time=2e21/flux_He,
)
F += -source * test_func[0][1] * dx
for m in range(2):
    for n in range(N + 1):
        if (n, m) != (0, 0):
            F += (sols[m][n] - prev_sols[m][n]) / dt * test_func[m][n] * dx + diff[m][
                n
            ] * dot(grad(sols[m][n]), grad(test_func[m][n])) * dx

for reac in reactions:
    prod = 1
    product_m = 0
    product_n = 0
    for s in reac["reactives"]:
        prod *= sols[s[0]][s[1]]
        product_m += s[0]
        product_n += s[1]

    for s in reac["reactives"]:
        F += -(-reac["k+"] * prod * test_func[s[0]][s[1]] * dx)
        F += -(reac["k-"] * sols[product_m][product_n] * test_func[s[0]][s[1]] * dx)

    F += -(
        -reac["k-"] * sols[product_m][product_n] * test_func[product_m][product_n] * dx
    )
    F += -(reac["k+"] * prod * test_func[product_m][product_n] * dx)


du = TrialFunction(c.function_space())
J = derivative(F, c, du)  # Define the Jacobian
set_log_level(30)

n = FacetNormal(mesh)
fig_num = 2

compute = True
ramp = 14
fluence = 2e21
flux_He = fluxflux # laflux
t_imp = fluence / flux_He
t_rest = 25
c_n.assign(c_n_back)
compute = True
ret = []
time = []
temp = []
flux = []
ret = []
t = 0
# density.assign(d)  # can't remember what's that for but d is not defined anywhere
T.assign(300)
source.implantation_time = t_imp
source.flux_He = flux_He
print("Total time of simulation: " + str(t_imp + t_rest + 2200 / ramp) + " s")
print(t_imp,t_rest,2200/ramp)
while t < t_imp + t_rest + 2200 / ramp:
    print(t, end="\r")
    t += float(dt)
    source.t = t
    if t > t_imp + t_rest:
        T.assign(300 + ramp * (t - (t_imp + t_rest)))
        dt.assign(0.5)

    problem = NonlinearVariationalProblem(F, c, bcs, J)
    solver = NonlinearVariationalSolver(problem)
    solver.parameters["newton_solver"]["absolute_tolerance"] = 1e9
    solver.parameters["newton_solver"]["relative_tolerance"] = 1e-10
    solver.parameters["newton_solver"]["maximum_iterations"] = 16
    nb_it, converged = solver.solve()
    if t < t_imp + t_rest:
        solving.adaptive_stepsize(nb_it, converged, dt, 1e-8, 1.05, t)
    else:
        solving.adaptive_stepsize(nb_it, converged, dt, 1e-8, 1.00, t)
    if t < t_imp and t + float(dt) > t_imp + 1:
        #print("on cut")
        #dt.assign(t_imp - t)
        dt.assign(1.0)
        print("tcut  "+str(t+float(dt))+"  s")
    c_n.assign(c)
    res = list(c.split())
    if t < t_imp + 2 and t > t_imp - 2 and compute is True:
        ret.append(assemble(res[N] * dx))
        ret.append(assemble(res[N + 1] * dx))
        ret.append(assemble(res[N + 2] * dx))
        ret.append(assemble(res[N + 3] * dx))
        ret.append(assemble(res[N + 4] * dx))
        plt.figure(fig_num)
        plt.bar(["V1", "He1V1", "He2V1", "He3V1", "He4V1"], ret)
        plt.ylim(top=6e16)
        plt.ylabel("Inventory (m-3)")
        plt.savefig("ret" + "{:.2e}".format(fluence) + ".png")
        compute = False
        fig_num += 1
    for i in range(0, nb_clusters):
        res[i].rename(str(i + 1), str(i + 1))
        files[i].write(res[i], t)
    if t > t_imp + t_rest:
        val = -assemble(diff[0][1] * dot(grad(res[0]), n) * ds)
        val += -assemble(diff[0][2] * dot(grad(res[1]), n) * ds)
        val += -assemble(diff[0][3] * dot(grad(res[2]), n) * ds)
        flux.append(val)
        time.append(t)
        temp.append(float(T))
plt.figure(1)
plt.plot(temp, flux, label="{:.2e}".format(ramp) + " K/s")

plt.figure(1)
plt.xlabel("Temperature (K)")
plt.title(flux_He)
plt.ylim(top=0.8e17,bottom=0)
plt.legend()
plt.ylabel("Desorption flux (He/m2/s)")
plt.savefig("plot_kornelsen"+str(flux_He)+".png")
