import fenics as f
import sympy as sp
import numpy as np
import FESTIM
import csv
from os import path


def find_maximum_loc(u, maximum, dof_coord):
    index = np.where(u.vector() == maximum)[0]
    x = dof_coord[index[0]][0]
    return x


def main(
        mesh_parameters, temperature, source, dt, t_final, nb_clusters=6,
        folder="temp"):
    # Files
    files = []
    for i in range(0, nb_clusters):
        files.append(
            f.XDMFFile(folder + "/c_" + str(i + 1) + ".xdmf"))

    files.append(f.XDMFFile(folder + "/cb.xdmf"))
    files.append(f.XDMFFile(folder + "/i.xdmf"))
    files.append(f.XDMFFile(folder + "/retention.xdmf"))
    for file_ in files:
        file_.parameters["flush_output"] = True
        file_.parameters["rewrite_function_mesh"] = False

    # Defining mesh

    size = mesh_parameters["size"]
    mesh = FESTIM.meshing.mesh_and_refine(mesh_parameters)
    n = f.FacetNormal(mesh)
    vm, sm = FESTIM.meshing.subdomains(
        mesh,
        {
            "mesh_parameters": mesh_parameters,
            "materials": [{"borders": [0, size], "id": 1}]}
        )
    dx = f.Measure('dx', domain=mesh)

    dt = f.Constant(dt)
    # Defining functions
    V = f.VectorFunctionSpace(mesh, 'P', 1, nb_clusters+2)
    V_DG1 = f.FunctionSpace(mesh, "DG", 1)
    V_DG0 = f.FunctionSpace(mesh, "DG", 0)
    dof_coord = V.tabulate_dof_coordinates()
    c = f.Function(V)
    v = f.TestFunction(V)
    sols = list(f.split(c))
    test_func = list(f.split(v))

    ini = []
    # # Uncomment the following to add initial conditions
    # ini = [{"value": 1e15*1e6*(FESTIM.x < 1)}]
    ini = [{"value": nb_clusters+1, "component": nb_clusters+1}]
    c_n = \
        FESTIM.initialising.initialise_solutions({"initial_conditions": ini}, V)
    prev_sols = f.split(c_n)
    bcs = []
    for i in range(0, nb_clusters+1):
        bcs.append(f.DirichletBC(V.sub(i), f.Constant(0), sm, 1))
        bcs.append(f.DirichletBC(V.sub(i), f.Constant(0), sm, 2))

    # Defining form
    if isinstance(temperature, (int, float)):
        T = f.Constant(temperature)
    else:
        T = f.Expression(
            sp.printing.ccode(temperature), t=0, degree=2)
    immobile_cluster_threshold = 7
    diff = [0 for i in range(nb_clusters + 1)]
    diff[0] = 2.95e-8*f.exp(-0.13/FESTIM.k_B/T)
    diff[1] = 3.24e-8*f.exp(-0.2/FESTIM.k_B/T)
    diff[2] = 2.26e-8*f.exp(-0.25/FESTIM.k_B/T)
    diff[3] = 1.68e-8*f.exp(-0.2/FESTIM.k_B/T)
    diff[4] = 0.520e-8*f.exp(-0.12/FESTIM.k_B/T)
    diff[5] = 0.120e-8*f.exp(-0.3/FESTIM.k_B/T)

    R1 = 3e-10
    R = []
    D = []

    def radius_pure_helium(i):
        a_0 = 0.318e-9
        pi = np.pi
        r = R1
        if i > 1:
            r += (3/(4*pi)*a_0**3/10*i)**(1/3)
            r -= (3/(4*pi)*a_0**3/10)**(1/3)
        return r
    # # Standard model
    atomic_density = 6.3e28
    # this has to be coherent with nb_clusters
    # C Bequart's values
    dissociation_energies = [
        None,  # He -> He + He  dum
        1,  # He2 -> He + He
        1.5,  # He3 -> He2 + He
        1.5,
        1.6,
        2,
        None,  # cb
    ]

    # S Blondel's values
    dissociation_energies = [
        None,
        0.8600000000000012,
        1.2399999999999984,
        1.5,
        1.0499999999999972,
        2.0100000000000016,
        None]
    reactions = []
    for i in range(0, nb_clusters):
        k_plus = 4*f.pi*(R1 + radius_pure_helium(i+1))*(diff[0] + diff[i])

        if dissociation_energies[i+1] is not None:
            k_minus = \
                atomic_density*k_plus * \
                f.exp(-dissociation_energies[i+1]/FESTIM.k_B/T)
        else:
            k_minus = 0
        reaction = {
            "reactives": [0, i],
            "k+": k_plus,
            "k-": k_minus,
        }
        reactions.append(reaction)

    R = [0]*(nb_clusters+1)
    D = [0]*(nb_clusters+1)
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
    if not isinstance(source, f.Function):
        source_expr = f.Expression(sp.printing.ccode(source), degree=1, t=0)
    else:
        source_expr = source
    F += -source_expr*test_func[0]*dx

    # # Extended model
    cb = sols[-2]
    cb_n = prev_sols[-2]
    av_i = sols[-1]
    av_i_n = prev_sols[-1]

    def radius(i):
        a_0 = 0.318e-9
        pi = np.pi
        return ((3**0.5)/4*a_0) + \
            (3/(4*pi) * (a_0**3)/2 * abs(i)/4)**(1/3) - \
            (3/(4*pi) * (a_0**3)/2)**(1/3)

    rb = radius(av_i_n)
    k_b_av = 4*f.pi*diff[0]*(R1 + rb)
    R[0] += - k_b_av*sols[0]*cb

    for i in range(0, nb_clusters + 1):
        # print(R[i])
        F += (sols[i] - prev_sols[i])/dt*test_func[i]*dx + \
            diff[i]*f.dot(f.grad(sols[i]), f.grad(test_func[i]))*dx
        F += (- R[i] - D[i])*test_func[i]*dx

    # doesn't work: d(cb*ib)/dt
    # F += (av_i*cb-av_i_n*cb_n)/dt *\
    #     test_func[-1]*dx

    # d(cb*ib)/dt = cb*di/dt + i*dcb/dt
    F += av_i*(cb-cb_n)/dt *\
        test_func[-1]*dx
    F += (cb+1e-5)*(av_i - av_i_n)/dt *\
        test_func[-1]*dx  # cb + 1e-5 is needed cause crashes if zero

    # reaction terms for i_b
    F += -((nb_clusters + 1)*R[-2])*test_func[-1]*dx
    F += -(k_b_av*sols[0]*cb)*test_func[-1]*dx

    du = f.TrialFunction(c.function_space())
    J = f.derivative(F, c, du)  # Define the Jacobian

    # Solving
    t = 0
    f.set_log_level(30)

    derived_quantities_global = [
        [
            't(s)', 'inventory(He/m2)', 'flux_surface_left(He/m2/s)',
            'max_ib(He)', 'x_max_ib(m)', 'mean_ib(He)', 'total_bubbles']]
    while t < t_final:
        print(format(t, "5.2e") + " s", end="\r")
        t += float(dt)
        source_expr.t = t
        T.t = t
        # solve
        problem = f.NonlinearVariationalProblem(F, c, bcs, J)
        solver = f.NonlinearVariationalSolver(problem)
        solver.parameters["newton_solver"]["absolute_tolerance"] = 1e10
        solver.parameters["newton_solver"]["relative_tolerance"] = 1e-10
        solver.parameters["newton_solver"]["maximum_iterations"] = 30
        nb_it, converged = solver.solve()

        # Post processing
        res = list(c.split())
        ret = 0
        for i in range(0, len(res)):
            if i == len(res) - 2:
                name = "cb"
            elif i == len(res) - 1:
                name = "ib"
            else:
                name = str(i+1)
            res[i].rename(name, name)
            files[i].write(res[i], t)
        immobile_clusters = sum(
            [(i+1)*res[i] for i in
                range(immobile_cluster_threshold-1, len(res)-2)]
            )
        retention = f.project(
            res[-1]*res[-2] + immobile_clusters)
        retention.rename('retention', 'retention')
        files[-1].write(retention, t)

        flux_left = 0
        for i in range(0, nb_clusters):
            if i < immobile_cluster_threshold:
                flux_left += f.assemble(diff[i]*f.dot(f.grad(res[i]), n)*f.ds)
        ib_max = \
            FESTIM.post_processing.calculate_maximum_volume(res[-1], vm, 1)
        x_max_ib = find_maximum_loc(res[-1], ib_max, dof_coord)
        immobile_clusters = \
            [res[i] for i in range(immobile_cluster_threshold-1, len(res)-1)]
        total_bubbles = f.assemble((sum(immobile_clusters))*dx)
        if total_bubbles > 0:
            mean_ib = f.assemble(retention*dx)/total_bubbles
        else:
            mean_ib = nb_clusters + 1
        derived_quantities_global.append(
            [t, f.assemble(retention*dx), flux_left, ib_max, x_max_ib, mean_ib,
             total_bubbles])

        # update
        FESTIM.solving.adaptive_stepsize(
            nb_it, converged, dt, 1e-8, 1.1, t)
        c_n.assign(c)
    FESTIM.export.write_to_csv(
        {"file": "derived_quantities.csv", "folder": folder},
        derived_quantities_global)


if __name__ == "__main__":
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
    center = 1.5e-9
    width = 1e-9
    distribution = 1/(width*(2*3.14)**0.5) * \
        sp.exp(-0.5*((FESTIM.x-center)/width)**2)
    flux = 1e22  # flux in He m-2 s-1
    source = distribution*flux/0.93
    main(
        mesh_parameters, dt=0.000001, t_final=50, temperature=1000,
        source=source)
