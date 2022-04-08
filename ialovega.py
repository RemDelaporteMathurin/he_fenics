from fenics import *
import sympy as sp
import FESTIM
import numpy as np
import matplotlib.pyplot as plt


def run_sim(source, T, t_final, t_stop):
    N = 4
    nb_clusters = 2*N + 1
    # Files
    files = []
    folder = "ialovega/"
    for i in range(0, nb_clusters):
        if i <= N - 1:
            m = 0
            n = i + 1
        else:
            m = 1
            n = i - N
        files.append(XDMFFile(folder + "/c_" + str(n) + "He" + str(m) + "V.xdmf"))
        files[i].parameters["flush_output"] = True
        files[i].parameters["rewrite_function_mesh"] = False

    # Defining mesh
    size = 100e-5

    mesh_parameters = {
        "size": size,
        "initial_number_of_cells": 300,
        "refinements": [
            {
                "x": 2e-6,
                "cells": 500
            },
            {
                "x": 100-9,
                "cells": 2000
            },
        ]
    }
    mesh = FESTIM.meshing.mesh_and_refine(mesh_parameters)
    vm, sm = FESTIM.meshing.subdomains(
        mesh, {"mesh_parameters": mesh_parameters,
            "materials": [{"borders": [0, size], "id": 1}]})
    dt = Constant(2)
    dx = Measure('dx', domain=mesh)
    # Defining functions
    V = VectorFunctionSpace(mesh, 'P', 1, nb_clusters)
    V_DG1 = FunctionSpace(mesh, "DG", 1)
    V_DG0 = FunctionSpace(mesh, "DG", 0)
    c = Function(V)
    v = TestFunction(V)
    sols = list(split(c))
    sols.insert(0, 0)
    sols = np.reshape(sols, (-1, int(len(sols)/2)))
    test_func = list(split(v))
    test_func.insert(0, 0)
    test_func = np.reshape(test_func, (-1, int(len(test_func)/2)))
    ini = []

    # Initial conditions
    width = 1.27e-9
    center = 1.44e-9
    distribution = 1#1/(width*(2*3.14)**0.5) * sp.exp(-0.5*((FESTIM.x-center)/width)**2)
    prm = {"initial_conditions": [{"value": 1e20, "component": N}]}
    c_n = FESTIM.initialising.initialise_solutions(prm, V)
    c_n_back = Function(V)
    c_n_back.assign(c_n)
    prev_sols = list(split(c_n))
    prev_sols.insert(0, 0)
    prev_sols = np.reshape(prev_sols, (-1, int(len(prev_sols)/2)))
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
    diff = np.zeros((2, N + 1)).tolist()
    diff[0][1] = 2.95e-8*exp(-0.13/FESTIM.k_B/T)
    diff[0][2] = 3.24e-8*exp(-0.2/FESTIM.k_B/T)
    diff[0][3] = 2.26e-8*exp(-0.25/FESTIM.k_B/T)
    diff[0][4] = 1.68e-8*exp(-0.2/FESTIM.k_B/T)
    diff[1][0] = 177e-8*exp(-1.29/FESTIM.k_B/T)  # Diffusivity of V1

    # Capture radii
    a_0 = 0.318e-9
    r_He1 = 3e-10
    r_V1 = 3**0.5/4*a_0
    r = [[0, r_He1], [r_V1]]

    i = 1
    while len(r[0]) != N + 1:
        i += 1
        r[0].append(r_He1 + (3/4/pi * a_0**3/10*i)**(1/3)-(3/4/pi * a_0**3/10)**(1/3))
    i = 1
    while len(r[1]) != N + 1:
        i += 1
        r[1].append(r_V1 + (3/4/pi * a_0**3/2*1)**(1/3)-(3/4/pi * a_0**3/2)**(1/3))


    # Reaction rates

    # Dissociation energies (emission of 1 He)
    E_b_He = [[None, None, 1, 1.1, 1.2, 1.3],
        [None, 4, 3, 2.75, 2.45, 2.3]]

    density = Constant(6.3e28)
    reactions = [
        # pure He cluster
        {
            "reactives": [(0, 1), (0, 1)],
            "k+": 4*pi*(diff[0][1] + diff[0][1])*(r[0][1] + r[0][1]),
            "k-": 4*pi*(diff[0][1] + diff[0][1])*(r[0][1] + r[0][1])*density*exp(-E_b_He[0][2]/FESTIM.k_B/T),
        },
        {
            "reactives": [(0, 1), (0, 2)],
            "k+": 4*pi*(diff[0][1] + diff[0][2])*(r[0][1] + r[0][2]),
            "k-": 4*pi*(diff[0][1] + diff[0][2])*(r[0][1] + r[0][2])*density*exp(-E_b_He[0][3]/FESTIM.k_B/T),
        },
        {
            "reactives": [(0, 1), (0, 3)],
            "k+": 4*pi*(diff[0][1] + diff[0][3])*(r[0][1] + r[0][3]),
            "k-": 4*pi*(diff[0][1] + diff[0][3])*(r[0][1] + r[0][3])*density*exp(-E_b_He[0][4]/FESTIM.k_B/T),
        },
        # absorbtion of 1 He1
        {
            "reactives": [(0, 1), (1, 0)],
            "k+": 4*pi*(diff[0][1] + diff[1][0])*(r[0][1] + r[1][0]),
            "k-": 4*pi*(diff[0][1] + diff[1][0])*(r[0][1] + r[1][0])*density*exp(-E_b_He[1][1]/FESTIM.k_B/T),
        },
        {
            "reactives": [(0, 1), (1, 1)],
            "k+": 4*pi*(diff[0][1] + diff[1][1])*(r[0][1] + r[1][1]),
            "k-": 4*pi*(diff[0][1] + diff[1][1])*(r[0][1] + r[1][1])*density*exp(-E_b_He[1][2]/FESTIM.k_B/T),
        },
        {
            "reactives": [(0, 1), (1, 2)],
            "k+": 4*pi*(diff[0][1] + diff[1][2])*(r[0][1] + r[1][2]),
            "k-": 4*pi*(diff[0][1] + diff[1][2])*(r[0][1] + r[1][2])*density*exp(-E_b_He[1][3]/FESTIM.k_B/T),
        },
        {
            "reactives": [(0, 1), (1, 3)],
            "k+": 4*pi*(diff[0][1] + diff[1][3])*(r[0][1] + r[1][3]),
            "k-": 4*pi*(diff[0][1] + diff[1][3])*(r[0][1] + r[1][3])*density*exp(-E_b_He[1][4]/FESTIM.k_B/T),
        },
        # absorption of 1 V1
        {
            "reactives": [(0, 2), (1, 0)],
            "k+": 4*pi*(diff[0][2] + diff[1][0])*(r[0][1] + r[1][0]),
            "k-": 0,#4*pi*(diff[0][2] + diff[1][0])*(r[0][1] + r[1][0])*density*exp(-E_b[1][2]/FESTIM.k_B/T),
        },
        {
            "reactives": [(0, 3), (1, 0)],
            "k+": 4*pi*(diff[0][3] + diff[1][0])*(r[0][1] + r[1][0]),
            "k-": 0,#4*pi*(diff[0][3] + diff[1][0])*(r[0][1] + r[1][0])*density*exp(-E_b[1][3]/FESTIM.k_B/T),
        },
        {
            "reactives": [(0, 4), (1, 0)],
            "k+": 4*pi*(diff[0][4] + diff[1][0])*(r[0][1] + r[1][0]),
            "k-": 0,#4*pi*(diff[0][4] + diff[1][0])*(r[0][1] + r[1][0])*density*exp(-E_b[1][4]/FESTIM.k_B/T),
        },
    ]
    F = 0
    center = 2.57e-9
    width = 1.58e-9
    distribution = 1/(width*(2*3.14)**0.5) * \
        sp.exp(-0.5*((FESTIM.x-center)/width)**2)

    # source = Expression("t < implantation_time ?" + sp.printing.ccode(distribution) + "*flux_He: 0", degree=2, t=0,flux_He=flux, implantation_time=implantation_time)
    F += -source*test_func[0][1]*dx
    for m in range(2):
        for n in range(N+1):
            if (n, m) != (0, 0):
                F += (sols[m][n] - prev_sols[m][n])/dt*test_func[m][n]*dx + \
                    diff[m][n]*dot(grad(sols[m][n]), grad(test_func[m][n]))*dx

    for reac in reactions:
        prod = 1
        product_m = 0
        product_n = 0
        for s in reac["reactives"]:
            prod *= sols[s[0]][s[1]]
            product_m += s[0]
            product_n += s[1]

        for s in reac["reactives"]:
            F += -(-reac["k+"]*prod*test_func[s[0]][s[1]]*dx)
            F += -(reac["k-"]*sols[product_m][product_n]*test_func[s[0]][s[1]]*dx)

        F += -(-reac["k-"]*sols[product_m][product_n] *
            test_func[product_m][product_n]*dx)
        F += -(reac["k+"]*prod*test_func[product_m][product_n]*dx)


    du = TrialFunction(c.function_space())
    J = derivative(F, c, du)  # Define the Jacobian

    n = FacetNormal(mesh)
    fig_num = 2

    compute = True
    c_n.assign(c_n_back)
    compute = True
    ret = []
    time = []
    temp = []
    flux = []
    ret = []
    t = 0

    print("Total time of simulation: " + str(t_final) + " s")
    while t < t_final:
        print(t, end='\r')
        t += float(dt)
        T.t = t
        source.t = t

        problem = NonlinearVariationalProblem(F, c, bcs, J)
        solver = NonlinearVariationalSolver(problem)
        solver.parameters["newton_solver"]["absolute_tolerance"] = 1e10
        solver.parameters["newton_solver"]["relative_tolerance"] = 1e-10
        solver.parameters["newton_solver"]["maximum_iterations"] = 30
        nb_it, converged = solver.solve()
        FESTIM.solving.adaptive_stepsize(nb_it, converged, dt, 1e-8, 1.1, t)
        if t >= t_stop:
            if float(dt) > 10:
                dt.assign(10)

        c_n.assign(c)
        res = list(c.split())
        # if t < t_imp + 2 and t > t_imp - 2 and compute is True:
        #     ret.append(assemble(res[N]*dx))
        #     ret.append(assemble(res[N+1]*dx))
        #     ret.append(assemble(res[N+2]*dx))
        #     ret.append(assemble(res[N+3]*dx))
        #     ret.append(assemble(res[N+4]*dx))
        #     plt.figure(fig_num)
        #     plt.bar(["V1", "He1V1", "He2V1", "He3V1", "He4V1"], ret)
        #     plt.ylim(top=6e16)
        #     plt.ylabel("Inventory (m-3)")
        #     plt.savefig("ret" + "{:.2e}".format(fluence) + ".png")
        #     compute = False
        for i in range(0, nb_clusters):
            res[i].rename(str(i+1), str(i+1))
            files[i].write(res[i], t)

        val = -assemble(diff[0][1]*dot(grad(res[0]), n)*ds)
        val += -assemble(diff[0][2]*dot(grad(res[1]), n)*ds)
        val += -assemble(diff[0][3]*dot(grad(res[2]), n)*ds)
        flux.append(val)
        time.append(t)
        temp.append(T(size/2))
    return np.array(flux), np.array(time), np.array(temp)


if __name__ == "__main__":

    t = FESTIM.t

    # create source
    center = 15e-9
    width = 4e-9
    distribution = 1/(width*(2*3.14)**0.5) * \
        sp.exp(-0.5*((FESTIM.x-center)/width)**2)
    flux = 8e14
    fluence = 8e16
    source_expr = distribution * flux * (t < fluence/flux)
    source = Expression(sp.printing.ccode(source_expr), degree=2, t=0)

    # create temperature expression
    ramp = 1
    T_expr = 300
    t_begin = fluence/flux
    t_end = fluence/flux
    T_maxs = [1250, 1350, 1450, 1550, 1650]
    for i in range(5):
        T_max = T_maxs[i]
        t_max = (T_max - 300)/ramp

        t_end += t_max
        T_expr += ramp*(t - t_begin)*(t_begin < t)*(t < t_end)
        t_begin += t_max

    T = Expression(sp.printing.ccode(T_expr), degree=0, t=0)
    set_log_level(30)

    flux_out, time, temp = run_sim(source=source, T=T, t_final=t_end, t_stop=fluence/flux - 20)
    plt.figure()
    t_beg = fluence/flux
    t_end = fluence/flux
    for i in range(5):
        T_max = T_maxs[i]
        t_max = (T_max - 300)/ramp
        t_end += t_max
        indexes = np.where((t_beg < time) & (time < t_end))
        plt.plot(temp[indexes], flux_out[indexes])
        t_beg += t_max
    plt.xlabel("T (K)")
    plt.ylabel("Flux He/m2/s")
    plt.savefig("ialovega/tds.png")
    plt.figure()
    plt.plot(time, temp)
    plt.xlabel("Time (s)")
    plt.ylabel("T (K)")
    plt.savefig("ialovega/temperature.png")
