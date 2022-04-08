from fenics import *
import sympy as sp
import numpy as np
import FESTIM
import csv
from os import path

nb_clusters = 6
# Files
folder = "temp"
files = []
for i in range(0, nb_clusters):
    files.append(
        XDMFFile("clustering/" + folder + "/c_" + str(i + 1) + ".xdmf"))
files.append(XDMFFile("clustering/" + folder + "/cb.xdmf"))
files.append(XDMFFile("clustering/" + folder + "/i.xdmf"))
for f in files:
    f.parameters["flush_output"] = True
    f.parameters["rewrite_function_mesh"] = False

# Defining mesh
size = 100e-6

mesh_parameters = {
    "size": size,
    "initial_number_of_cells": 900,
    "refinements": [
        {
            "x": 2e-6,
            "cells": 600
        },
        {
            "x": 50-9,
            "cells": 50
        },
    ]
}
mesh = FESTIM.meshing.mesh_and_refine(mesh_parameters)
vm, sm = FESTIM.meshing.subdomains(
    mesh, {"mesh_parameters": mesh_parameters,
           "materials": [{"borders": [0, size], "id": 1}]})
dx = Measure('dx', domain=mesh)

sample_size = 6
points = np.c_[np.random.uniform(low=700, high=780, size=(sample_size,)),
               np.random.uniform(low=20, high=21, size=(sample_size,))]
for point in points:
    temperature = point[0]
    flux = 10**point[1]
    print(temperature, flux)
    if not path.exists("bubble_size/{:.3e}_{:.2e}.csv".format(temperature, flux)):
        dt = Constant(0.000001)
        # Defining functions
        V = VectorFunctionSpace(mesh, 'P', 1, nb_clusters+2)
        V_DG1 = FunctionSpace(mesh, "DG", 1)
        V_DG0 = FunctionSpace(mesh, "DG", 0)
        c = Function(V)
        v = TestFunction(V)
        sols = list(split(c))
        test_func = list(split(v))

        ini = []
        # # Uncomment the following to add initial conditions
        # ini = [{"value": 1e15*1e6*(FESTIM.x < 1)}]
        ini = [{"value": 7, "component": nb_clusters+2-1}]
        c_n, prev_sols = \
            FESTIM.initialise_solutions.initialising_solutions(V, ini)

        bcs = []
        for i in range(0, nb_clusters+2-1):
            bcs.append(DirichletBC(V.sub(i), Constant(0), sm, 1))
            bcs.append(DirichletBC(V.sub(i), Constant(0), sm, 2))

        # Defining form
        T = Constant(temperature)
        diff = [0 for i in range(nb_clusters + 1)]
        diff[0] = 2.95e-8*exp(-0.13/FESTIM.k_B/T)
        diff[1] = 3.24e-8*exp(-0.2/FESTIM.k_B/T)
        diff[2] = 2.26e-8*exp(-0.25/FESTIM.k_B/T)
        diff[3] = 1.68e-8*exp(-0.2/FESTIM.k_B/T)
        diff[4] = 0.520e-8*exp(-0.12/FESTIM.k_B/T)
        diff[5] = 0.120e-8*exp(-0.3/FESTIM.k_B/T)

        R1 = 3e-10
        R = []
        D = []

        # # Standard model
        reactions = []
        for i in range(0, nb_clusters):
            reaction = {
                "reactives": [0, i],
                "k+": 4*pi*(R1 + (i+1)**(1/3)*R1)*(diff[0] + diff[i]),
                "k-": 0,
            }
            reactions.append(reaction)

        R = [0]*(nb_clusters+2)
        D = [0]*(nb_clusters+2)
        for reaction in reactions:
            species = reaction["reactives"]
            prod = 1
            for s in species:
                prod *= sols[s]
            for s in species:
                R[s] += -reaction["k+"]*prod
                D[s] += reaction["k-"]*sols[sum(species)+1]
            R[sum(species)+1] += reaction["k+"]*prod
            D[sum(species)+1] += -reaction["k-"] * sols[sum(species)+1]

        F = 0
        center = 1.5e-9
        width = 1e-9
        distribution = 1/(width*(2*3.14)**0.5) * \
            sp.exp(-0.5*((FESTIM.x-center)/width)**2)
        # flux = 1e21  # flux in He m-2 s-1
        source = Expression(sp.printing.ccode(distribution*flux/0.93), degree=1)
        F += -source*test_func[0]*dx

        # # Extended model
        cb = sols[-2]
        av_i = sols[-1]
        av_i_n = prev_sols[-1]
        a_0 = 0.318e-9
        rb = (3**0.5/4*a_0) + (3/(4*pi) * a_0**3/2 *abs(av_i_n)/4)**(1/3) - (3/(4*pi) * a_0**3/2)**(1/3)
        k_b_av = 4*pi*diff[0]*(R1 + rb)
        R[0] += - k_b_av*sols[0]*cb

        for i in range(0, nb_clusters + 1):
            # print(R[i])
            F += (sols[i] - prev_sols[i])/dt*test_func[i]*dx + \
                diff[i]*dot(grad(sols[i]), grad(test_func[i]))*dx
            F += (- R[i] - D[i])*test_func[i]*dx

        F += (cb + 1)*(av_i - av_i_n)/dt *\
            test_func[-1]*dx
        F += -(k_b_av*sols[0]*cb)*test_func[-1]*dx
        F += -(R[len(R) - 2]*(nb_clusters + 1 - av_i))*test_func[-1]*dx

        du = TrialFunction(c.function_space())
        J = derivative(F, c, du)  # Define the Jacobian

        # Solving
        t = 0
        ret_old = 0
        ret = 0
        set_log_level(30)
        av_i_file = XDMFFile("i.xdmf")
        # while t < 50:


        def radius(i):
            pi = np.pi
            return (3**0.5/4*a_0) + (3/(4*pi) * a_0**3/2 *abs(i)/4)**(1/3) - (3/(4*pi) * a_0**3/2)**(1/3)


        output = []
        max_i = 6
        while radius(max_i) < 10e-9:
            if t > 1e6:
                break
            print(format(t, "5.2e") + " s", end="\r")
            t += float(dt)
            problem = NonlinearVariationalProblem(F, c, bcs, J)
            solver = NonlinearVariationalSolver(problem)
            solver.parameters["newton_solver"]["absolute_tolerance"] = 1e10
            solver.parameters["newton_solver"]["relative_tolerance"] = 1e-10
            solver.parameters["newton_solver"]["maximum_iterations"] = 16
            nb_it, converged = solver.solve()
            res = list(c.split())
            # ret = 0
            # for i in range(0, nb_clusters): # [0, nb_clusters-1, nb_clusters-2]: 
            #     # ret += assemble((i+1)*res[i]*dx)
            #     res[i].rename(str(i+1), str(i+1))
            #     files[i].write(res[i], t)
            FESTIM.solving.adaptive_stepsize(nb_it, converged, dt, 1e-8, 1.1, t, 1e16, 1)
            c_n.assign(c)
            # av_i.assign(project(sols[nb_clusters-1]/(sols[nb_clusters-2]), V_DG1))
            # av_i_file.write(av_i, t)
            # ret_old = ret
            max_i = project(res[-1], V_DG0).vector().max()
            output.append([t, radius(max_i)])
            # print(max_i, radius(max_i))

        FESTIM.export.write_to_csv({"file": "bubble_size/{:.3e}_{:.2e}.csv".format(temperature, flux)}, output)
