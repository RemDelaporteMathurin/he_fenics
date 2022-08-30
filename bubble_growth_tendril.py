from fenics import *
import sympy as sp
import numpy as np
from src import meshing, initialising, post_processing, solving


k_B = 8.617e-5  # eV/K
R_g = 8.314  # J/mol/K


def find_maximum_loc(u, maximum, dof_coord):
    index = np.where(u.vector() == maximum)[0]
    x = dof_coord[index[0]][0]
    return x


def main(
    mesh_parameters,
    temperature,
    source,
    dt,
    t_final,
    nb_clusters=6,
    folder="./",
    soret=False,
    bursting=False,
):
    # Files
    files = []
    for i in range(0, nb_clusters):
        files.append(XDMFFile(folder + "/c_" + str(i + 1) + ".xdmf"))

    files.append(XDMFFile(folder + "/cb.xdmf"))
    files.append(XDMFFile(folder + "/m.xdmf"))
    files.append(XDMFFile(folder + "/i.xdmf"))
    files.append(XDMFFile(folder + "/T.xdmf"))
    files.append(XDMFFile(folder + "/retention.xdmf"))
    for f in files:
        f.parameters["flush_output"] = True
        f.parameters["rewrite_function_mesh"] = False

    # Defining mesh

    size = mesh_parameters["size"]
    mesh = meshing.mesh_and_refine(mesh_parameters)
    n = FacetNormal(mesh)
    vm, sm = meshing.subdomains(
        mesh,
        {
            "mesh_parameters": mesh_parameters,
            "materials": [{"borders": [0, size], "id": 1}],
        },
    )
    dx = Measure("dx", domain=mesh)

    dt = Constant(dt)
    # Defining functions
    V = VectorFunctionSpace(mesh, "P", 1, nb_clusters + 3)
    W = FunctionSpace(mesh, "P", 1)
    V_DG1 = FunctionSpace(mesh, "DG", 1)
    V_DG0 = FunctionSpace(mesh, "DG", 0)
    dof_coord = V.tabulate_dof_coordinates()
    c = Function(V)
    v = TestFunction(V)
    sols = list(split(c))
    test_func = list(split(v))

    ini = []
    # # Uncomment the following to add initial conditions
    # ini = [{"value": 1e15*1e6*(x < 1)}]
    ini = [
        {
            "value": 1,
            "component": nb_clusters + 1,
        },  # initialise vacancies in bubbles to 1
        {
            "value": nb_clusters + 1,
            "component": nb_clusters + 2,
        },  # initialise helium content in bubbles to nb_clusters + 1
    ]
    c_n = initialising.initialise_solutions({"initial_conditions": ini}, V)
    prev_sols = split(c_n)
    bcs = []
    for i in range(0, nb_clusters + 1):
        bcs.append(DirichletBC(V.sub(i), Constant(0), sm, 1))
        bcs.append(DirichletBC(V.sub(i), Constant(0), sm, 2))

    # Defining form
    if isinstance(temperature, (int, float)):
        T = Constant(temperature)
    else:
        ccode = sp.printing.ccode(temperature)
        T = Expression(ccode, t=0, degree=2)
        if "t" not in ccode:
            T = interpolate(T, W)

    immobile_cluster_threshold = 7
    diff = [0 for i in range(nb_clusters + 1)]
    diff[0] = 2.9e-8*exp(-0.13/k_B/T)
    diff[1] = 3.2e-8 * exp(-0.2 / k_B / T)
    diff[2] = 2.3e-8 * exp(-0.25 / k_B / T)
    diff[3] = 1.7e-8 * exp(-0.2 / k_B / T)
    diff[4] = 5.0e-9 * exp(-0.12 / k_B / T)
    diff[5] = 1.0e-9 * exp(-0.3 / k_B / T)

    free_enthalpy = [1 for _ in range(nb_clusters + 1)]
    entropy = [1 for _ in range(nb_clusters + 1)]

    R1 = 3e-10
    R = []
    D = []

    def radius_pure_helium(i):
        a_0 = 0.318e-9
        pi = np.pi
        r = R1
        if i > 1:
            r += (3 / (4 * pi) * a_0**3 / 10 * i) ** (1 / 3)
            r -= (3 / (4 * pi) * a_0**3 / 10) ** (1 / 3)
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
        None,
    ]
    reactions = []
    for i in range(0, nb_clusters):
        k_plus = 4 * pi * (R1 + radius_pure_helium(i + 1)) * (diff[0] + diff[i])

        if dissociation_energies[i + 1] is not None:
            k_minus = (
                atomic_density * k_plus * exp(-dissociation_energies[i + 1] / k_B / T)
            )
        else:
            k_minus = 0
        reaction = {
            "reactives": [0, i],
            "k+": k_plus,
            "k-": k_minus,
        }
        reactions.append(reaction)

    R = [0] * (nb_clusters + 1)
    D = [0] * (nb_clusters + 1)
    for reaction in reactions:
        species = reaction["reactives"]
        prod = 1
        for s in species:
            prod *= sols[s]
        for s in species:
            R[s] += -reaction["k+"] * prod
            D[s] += reaction["k-"] * sols[sum(species) + 1]
        R[sum(species) + 1] += reaction["k+"] * prod
        D[sum(species) + 1] += -reaction["k-"] * sols[sum(species) + 1]
    F = 0
    if not isinstance(source, Function):
        source_expr = Expression(sp.printing.ccode(source), degree=1, t=0)
    else:
        source_expr = source
    F += -source_expr * test_func[0] * dx

    # # Extended model
    cb = sols[-3]
    cb_n = prev_sols[-3]
    av_m = sols[-2]
    av_m_n = prev_sols[-2]
    av_i = sols[-1]
    av_i_n = prev_sols[-1]

    def radius(nb_V):
        """Returns the radius of a cluster based on the number of vacancies

        Args:
            nb_V (float, fenics.Function): the number of vacancies

        Returns:
            float or ufl.product: the radius in m
        """    
        a_0 = 0.318e-9
        pi = np.pi
        return (
            ((3**0.5) / 4 * a_0)
            + (3 / (4 * pi) * (a_0**3) / 2 * nb_V ) ** (1 / 3)   ###
            - (3 / (4 * pi) * (a_0**3) / 2) ** (1 / 3)
        )

    rb = radius(abs(av_m_n))
    k_b_av = 4 * pi * diff[0] * (R1 + rb)
    R[0] += -k_b_av * sols[0] * cb

    for i in range(0, nb_clusters + 1):
        # print(R[i])
        F += (sols[i] - prev_sols[i]) / dt * test_func[i] * dx + diff[i] * dot(
            grad(sols[i]), grad(test_func[i])
        ) * dx
        # soret
        if isinstance(T, Function) and soret:
            print("Adding soret effect")
            Q = free_enthalpy[i] * T + entropy[i]
            F += (
                dot(
                    #diff[i] * Q * sols[i] / (R_g * T**2) * grad(T),
                    diff[i] * -0.0065 * sols[i] * grad(T),
                    grad(test_func[i]),
                )
                * dx
            )
        F += (-R[i] - D[i]) * test_func[i] * dx

    # doesn't work: d(cb*ib)/dt
    # F += (av_i*cb-av_i_n*cb_n)/dt *\
    #     test_func[-1]*dx

    # d(cb*ib)/dt = cb*di/dt + i*dcb/dt
    F += av_i * (cb - cb_n) / dt * test_func[-1] * dx
    F += (
        (cb + 1e-6) * (av_i - av_i_n) / dt * test_func[-1] * dx
    )  # cb + 1e-5 is needed cause crashes if zero

    # reaction terms for i_b
    F += -((nb_clusters + 1) * R[-2]) * test_func[-1] * dx
    F += -(k_b_av * sols[0] * cb) * test_func[-1] * dx

    # bursting term for i_b
    f = Constant(2e3)
    # k_burst = alpha * Expression("exp(-x[0])", degree=2)
    x = SpatialCoordinate(mesh)[0]
    k_burst =  f * (1 - (x - rb) / x)
    if bursting:
        F += av_i * k_burst * cb * test_func[-1] * dx

     # d(cb*m)/dt = cb * dm/dt + m*dcb/dt
    F += av_m * (cb - cb_n) / dt * test_func[-2] * dx
    F += (
        (cb + 1e-6) * (av_m - av_m_n) / dt * test_func[-2] * dx
    )  # cb + 1e-5 is needed cause crashes if zero

    # reaction terms for m
    def max_nb_he(nb_vacancies):
        return nb_vacancies * 4

    X = (av_i_n + 1) / max_nb_he(av_m_n) - 1
    sigmoid = 1 / (
        1 + exp(-50 * X)
    )  # 50 is chosen to have a high derivative in the sigmoid
    F += -R[-2] * test_func[-2] * dx
    F += -sigmoid * (k_b_av * sols[0] * cb) * test_func[-2] * dx

    du = TrialFunction(c.function_space())
    J = derivative(F, c, du)  # Define the Jacobian

    # Solving
    t = 0
    set_log_level(30)

    derived_quantities_global = [
        [
            "t(s)",
            "inventory(He/m2)",
            "flux_surface_left(He/m2/s)",
            "max_ib(He)",
            "max_m(V)",
            "x_max_ib(m)",
            "mean_ib(He)",
            "mean_m(V)",
            "mean_ratio_HeV",
            "total_bubbles",
        ]
    ]
    while t < t_final:
        print(format(t, "5.2e") + " s", end="\r")
        t += float(dt)
        source_expr.t = t
        T.t = t
        # solve
        problem = NonlinearVariationalProblem(F, c, bcs, J)
        solver = NonlinearVariationalSolver(problem)
        solver.parameters["newton_solver"]["absolute_tolerance"] = 1e10
        solver.parameters["newton_solver"]["relative_tolerance"] = 1e-11
        solver.parameters["newton_solver"]["maximum_iterations"] = 30
        nb_it, converged = solver.solve()

        # Post processing
        res = list(c.split())
        for i in range(0, len(res)):
            if i == len(res) - 3:
                name = "cb"
            elif i == len(res) - 2:
                name = "m"
            elif i == len(res) - 1:
                name = "ib"
            else:
                name = str(i + 1)
            res[i].rename(name, name)
            #files[i].write(res[i], t)
        immobile_clusters = sum(
            [
                (i + 1) * res[i]
                for i in range(immobile_cluster_threshold - 1, len(res) - 3)
            ]
        )
        retention = project(res[-1] * res[-3] + immobile_clusters)
        retention.rename("retention", "retention")
        total_vacancies = project(res[-2] * res[-3])
        #files[-1].write(retention, t)
        #files[-2].write(T, t)
        flux_left = 0
        for i in range(0, nb_clusters):
            if i < immobile_cluster_threshold:
                flux_left += assemble(diff[i] * dot(grad(res[i]), n) * ds)
        ib_max = post_processing.calculate_maximum_volume(res[-1], vm, 1)
        m_max = post_processing.calculate_maximum_volume(res[-2], vm, 1)
        x_max_ib = find_maximum_loc(res[-1], ib_max, dof_coord)
        immobile_clusters = [
            res[i] for i in range(immobile_cluster_threshold - 1, len(res) - 2)
        ]
        total_bubbles = assemble((sum(immobile_clusters)) * dx)
        if total_bubbles > 0:
            mean_ib = assemble(retention * dx) / total_bubbles
            mean_m = assemble(total_vacancies* dx) / total_bubbles
            mean_ratio_he_v = mean_ib/mean_m  # TODO check that this is the real expression
        else:
            mean_ib = nb_clusters + 1
            mean_m = 1
        derived_quantities_global.append(
            [
                t,
                assemble(retention * dx),
                flux_left,
                ib_max,
                m_max,
                x_max_ib,
                mean_ib,
                mean_m,
                mean_ratio_he_v,
                total_bubbles,
            ]
        )

        # update
        solving.adaptive_stepsize(nb_it, converged, dt, 1e-8, 1.05, t)
        c_n.assign(c)

        if t + float(dt) > t_final:
            dt.assign(t_final - t)

    post_processing.write_to_csv(
        {"file": folder + "r-derived_quantities.csv"}, derived_quantities_global
    )
    # {"file": "r-derived_quantities.csv", "folder": folder},
    # derived_quantities_global)

    for i in range(len(res) - 2):
        post_processing.export_txt(folder + "r-{}".format(i + 1), res[i], W)
        # post_processing.export_txt(folder + "/r-{}".format(i+1), res[i], W)
    post_processing.export_txt(folder + "r-cb", res[-3], W)
    post_processing.export_txt(folder + "r-m", res[-2], W)
    post_processing.export_txt(folder + "r-ib", res[-1], W)
    if not isinstance(T, Function):
        T = project(T, W)
    post_processing.export_txt(folder + "r-T", T, W)
    post_processing.export_txt(folder + "radius", radius(abs(res[-2])), W)  # Rayon bulle <rb> papier


if __name__ == "__main__":
    x = sp.Symbol("x[0]")
    t = sp.Symbol("t")
    size = 30e-9

    mesh_parameters = {
        "size": size,
        "initial_number_of_cells":900,
        #"refinements": [
        #    {"x": 2e-6, "cells": 600},
        #    {"x": 50e-9, "cells": 100},
        #],
    }

    #   TS Xolotl W100
    distribution = (x <= 10e-9) *  1e9/16.39 *  (7.00876507 + 0.6052078 * x*1e9 - 3.01711048 *(x*1e9)**2 + 1.36595786 * (x*1e9)**3 - 0.295595 * (x*1e9)**4 + 0.03597462 * (x*1e9)**5 - 0.0025142 * (x*1e9)**6 + 0.0000942235 * (x*1e9)**7 - 0.0000014679 * (x*1e9)**8)
    flux = 1e22  # flux in He m-2 s-1


    # center = 1.5e-9
    # width = 1e-9
    # distribution = 1/(width*(2*3.14)**0.5) * sp.exp(-0.5*((x-center)/width)**2)
    # print('Integration of distribution:', float(sp.integrate(distribution, (x, 0.01, 30e-9))))
    # distribution = distribution/sp.integrate(distribution, (x, 0, 30e-9))

    source = distribution*flux

    main(
        mesh_parameters,
        dt=1.0,
        t_final=5011,
        temperature=500,
        folder="./",
        source=source,
        soret=False,
        bursting=True,
    )