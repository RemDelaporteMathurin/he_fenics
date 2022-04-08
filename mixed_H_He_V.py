from collections import Counter
import math
import sympy as sp
import csv

from fenics import *
import FESTIM

import sys
sys.setrecursionlimit(1500)

diffusion_coefficients = [
    # nb_Vac nb_H nb_He

    # 0 vac
    [
        # 0 H
        [[None, None], [2.95e-8, 0.13], [3.24e-8, 0.2], [2.26e-8, 0.25], [1.68e-8, 0.2], [0.52e-8, 0.12]],
        # 1 H
        [[1.9e-7, 0.2], [2.95e-8, 0.13], [3.24e-8, 0.2], [2.26e-8, 0.25], [1.68e-8, 0.2], [0.52e-8, 0.12]],
        # 2 H
        [[1.9e-7, 0.2], [2.95e-8, 0.13], [3.24e-8, 0.2], [2.26e-8, 0.25], [1.68e-8, 0.2], [0.52e-8, 0.12]],
        # 3 H
        [[1.9e-7, 0.2], [2.95e-8, 0.13], [3.24e-8, 0.2], [2.26e-8, 0.25], [1.68e-8, 0.2], [0.52e-8, 0.12]],
        # 4 H
        [[1.9e-7, 0.2], [2.95e-8, 0.13], [3.24e-8, 0.2], [2.26e-8, 0.25], [1.68e-8, 0.2], [0.52e-8, 0.12]],
        # 5 H
        [[1.9e-7, 0.2], [2.95e-8, 0.13], [3.24e-8, 0.2], [2.26e-8, 0.25], [1.68e-8, 0.2], [0.52e-8, 0.12]],
    ],
    # 1 vac
    [
        # 0 H
        [0],
        # 1 H
        [0]
    ]

]


def find_cluster(clusters, n_V, n_H, n_He):
    for cluster in clusters:
        if (cluster.n_V, cluster.n_H, cluster.n_He) == (n_V, n_H, n_He):
            return cluster
    print("Couldn't find cluster {}V.{}H.{}He".format(n_V, n_H, n_He))


class Cluster():
    def __init__(self, n_V, n_H, n_He):
        self.n_He = n_He
        self.n_H = n_H
        self.n_V = n_V
        self.D_0, self.E_D = self.compute_diffusion_coeff(n_V, n_H, n_He)
        self.E_b_H, self.E_b_He, self.E_b_V = \
            self.compute_binding_energies(n_V, n_H, n_He)
        self.reactions = self.compute_reactions(n_V, n_H, n_He)
        self.concentration = None
        self.previous_concentration = None
        self.test_function = None
        self.radius = self.compute_radius(n_V, n_H, n_He)
        self.index = None
        self.name = "{}V.{}H.{}He".format(n_V, n_H, n_He)

    def get_index(self):
        return self.index

    def compute_diffusion_coeff(self, n_V, n_H, n_He):
        if n_V != 0:
            return 0, 0
        else:
            try:
                D = diffusion_coefficients[n_V][n_H][n_He]
                D_0, E_D = D[0], D[1]
            except IndexError:
                print(
                    'Diffusion coefficient not found in table. 0 m2/s assumed.'
                    )
                D_0, E_D = 0, 0
                pass
            return D_0, E_D

    def compute_radius(self, n_V, n_H, n_He):
        a_0 = 0.318e-9
        rHe1 = 0.3e-9
        rHe0V1 = 3**0.5/4*a_0
        if n_V == 0:
            radius = rHe1 + (3/(4*math.pi)*a_0**3/10 * n_He)**(1/3) - \
                (3/(4*math.pi)*a_0**3/10)**(1/3)
        else:
            radius = rHe0V1 + (3/(4*math.pi)*a_0**3/2 * n_V)**(1/3) - \
                (3/(4*math.pi)*a_0**3/2)**(1/3)
        return radius

    def compute_binding_energies(self, n_V, n_H, n_He):
        vals = []
        with open('binding_energies.csv') as csvfile:
            reader = csv.reader(csvfile, delimiter=";")
            next(reader)
            for row in reader:
                n_V_ = int(row[0])
                n_He_ = int(row[1])
                n_H_ = int(row[2])
                if (n_V_, n_He_, n_H_) == (n_V, n_He, n_H):
                    for i in range(3):
                        if row[3 + i] == "None":
                            vals.append(None)
                        elif row[3 + i] == "":
                            vals.append(0.2)
                        else:
                            vals.append(float(row[3 + i]))
                    break
        energy_H = vals[0]
        energy_He = vals[1]
        energy_V = vals[2]
        return energy_H, energy_He, energy_V

    def compute_reactions(self, n_V, n_H, n_He):
        reactants = []
        if (n_V, n_H, n_He) not in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:

            if n_V != 0:
                if not(n_H > 1 and n_He == 0):  # ex: reaction V + H6 --> VH6 not possible
                    reactants.append([])
                    reactants[-1].append((n_V - 1, n_H, n_He))
                    reactants[-1].append((1, 0, 0))

            if n_H != 0:
                reactants.append([])
                reactants[-1].append((n_V, n_H - 1, n_He))
                reactants[-1].append((0, 1, 0))
            if n_He != 0:
                if not(n_H > 1 and n_V == 0):  # ex: reaction He + H2 -> HeH2 not possible
                    reactants.append([])
                    reactants[-1].append((n_V, n_H, n_He - 1))
                    reactants[-1].append((0, 0, 1))

        # delete copies
        ctr = Counter(frozenset(x) for x in reactants)
        b = [ctr[frozenset(x)] == 1 for x in reactants]
        unique_reactants = []
        for i in range(len(reactants)):
            if b[i] is True:
                unique_reactants.append(reactants[i])
            else:
                unique_reactants.append(reactants[0])
                break
        return unique_reactants


class Simulation():
    def __init__(self, n_He_max, n_H_max, n_Vac_max, parameters={}):
        self.n_He_max, self.n_H_max, self.n_Vac_max = \
            n_He_max, n_H_max, n_Vac_max
        self.parameters = parameters
        self.clusters = self.create_clusters()
        self.mesh = None
        self.t = 0
        self.dt = Constant(parameters["solving"]["initial_stepsize"])
        self.T = parameters["temperature"]
        self.density = 6.3e28
        self.t_final = parameters["solving"]["final_time"]
        self.F = 0
        self.u, self.u_n, self.v = None, None, None
        self.bcs = []
        self.file_out = XDMFFile("Solution/out.xdmf")
        self.file_H_tot = XDMFFile("Solution/H_tot.xdmf")
        self.file_He_tot = XDMFFile("Solution/He_tot.xdmf")
        self.file_V_tot = XDMFFile("Solution/V_tot.xdmf")
        self.file_VHe_tot = XDMFFile("Solution/VHe_tot.xdmf")
        self.file_VH_tot = XDMFFile("Solution/VH_tot.xdmf")
        self.file_VHeH_tot = XDMFFile("Solution/VHeH_tot.xdmf")
        self.derived_quantities = []

    def create_clusters(self):
        clusters = []

        for n_V in range(self.n_Vac_max + 1):
            for n_H in range(self.n_H_max + 1):
                for n_He in range(self.n_He_max + 1):
                    if (n_V, n_H, n_He) != (0, 0, 0) and n_H + n_He <= 6:
                        if not (n_V == 0 and n_He == 0 and n_H > 1):
                            cluster = Cluster(n_V, n_H, n_He)
                            clusters.append(cluster)
        return clusters

    def attribute_functions_to_clusters(self):
        solutions = list(split(self.u))
        prev_solutions = list(split(self.u_n))
        test_functions = list(split(self.v))
        i = 0
        for cluster in self.clusters:
            cluster.concentration = solutions[i]
            cluster.previous_concentration = prev_solutions[i]
            cluster.test_function = test_functions[i]
            cluster.index = i
            i += 1

    def create_formulation(self):
        k_B = FESTIM.k_B
        F = 0
        clusters = self.clusters
        dt = self.dt
        T = self.T
        density = self.density
        for cluster in clusters:
            if (cluster.n_V, cluster.n_H, cluster.n_He) != (0, 0, 0):
                F += -(
                    cluster.concentration - cluster.previous_concentration) / \
                        dt*cluster.test_function*dx
                if cluster.D_0 != 0:
                    F += -cluster.D_0*exp(-cluster.E_D/(k_B*T)) * \
                        dot(
                            grad(cluster.concentration),
                            grad(cluster.test_function))*dx
            for reac in cluster.reactions:
                r1 = find_cluster(clusters, reac[0][0], reac[0][1], reac[0][2])
                r2 = find_cluster(clusters, reac[1][0], reac[1][1], reac[1][2])
                D1, D2 = r1.D_0 * exp(-r1.E_D/(k_B*T)), \
                    r2.D_0 * exp(-r2.E_D/(k_B*T))
                k_plus = 4*pi*(r1.radius + r1.radius)*(D1 + D2)
                if (r2.n_V, r2.n_H, r2.n_He) == (0, 1, 0):
                    binding_energy = cluster.E_b_H
                elif (r2.n_V, r2.n_H, r2.n_He) == (0, 0, 1):
                    binding_energy = cluster.E_b_He
                elif (r2.n_V, r2.n_H, r2.n_He) == (1, 0, 0):
                    binding_energy = cluster.E_b_V
                else:
                    raise ValueError('Reactant 2 must be He, H or V')
                # reaction
                F += k_plus*r1.concentration * \
                    r2.concentration*cluster.test_function*dx
                F -= k_plus*r1.concentration * \
                    r2.concentration*r1.test_function*dx
                F -= k_plus*r1.concentration * \
                    r2.concentration*r2.test_function*dx
                # dissociation
                k_minus = k_plus * density * exp(-binding_energy/(k_B*T))
                F -= k_minus*cluster.concentration*cluster.test_function*dx
                F += k_minus*cluster.concentration*r1.test_function*dx
                F += k_minus*cluster.concentration*r2.test_function*dx
        return F

    def create_bcs(self, V, sm):
        bcs = []
        for i in range(0, V.num_sub_spaces()):
            bcs.append(DirichletBC(V.sub(i), Constant(0), sm, 1))
            bcs.append(DirichletBC(V.sub(i), Constant(0), sm, 2))
        return bcs

    def create_source_terms(self):
        source_form = 0
        for source in self.parameters["sources"]:
            source_val = Expression(
                sp.printing.ccode(source["distribution"]*source["value"]),
                degree=2, t=0)
            n_V, n_H, n_He = \
                source["cluster"]["n_V"], \
                source["cluster"]["n_H"], \
                source["cluster"]["n_He"]
            cluster = find_cluster(self.clusters, n_V=n_V, n_H=n_H, n_He=n_He)
            source_form += source_val*cluster.test_function*dx

        return source_form

    def run(self):

        # Defining mesh
        self.mesh = FESTIM.meshing.mesh_and_refine(self.parameters["mesh"])
        vm, sm = FESTIM.meshing.subdomains(
            self.mesh, {"mesh_parameters": mesh_parameters,
                        "materials": [
                            {"borders": [0, self.parameters["mesh"]["size"]],
                             "id": 1}]})
        dx = Measure('dx', domain=self.mesh)
        V = VectorFunctionSpace(
            self.mesh, 'P', 1,
            len(self.clusters))
        W = FunctionSpace(self.mesh, 'P', 1)
        self.u, self.u_n = Function(V), Function(V)

        self.v = TestFunction(V)

        self.attribute_functions_to_clusters()
        ini = []
        for initial_conc in [(1, 0, 0)]:
            cluster = find_cluster(self.clusters, *initial_conc)
            ini.append({"value": 2e25, "component": cluster.get_index()})
        self.u_n.assign(
            FESTIM.initialise_solutions.initialising_solutions(V, ini)[0])

        self.bcs = self.create_bcs(V, sm)

        self.F = self.create_formulation()
        self.F += self.create_source_terms()
        self.time_stepping()
        FESTIM.export.write_to_csv(
            {"file": "Solution/out.csv"}, self.derived_quantities)

    def time_stepping(self):
        # W = FunctionSpace(self.mesh, 'P', 1)
        du = TrialFunction(self.u.function_space())
        J = derivative(self.F, self.u, du)  # Define the Jacobian
        set_log_level(30)
        self.derived_quantities.append(["t", "T", "inv_H", "inv_He"])
        while self.t < self.t_final:
            print(format(self.t, "5.2e") + " s", end="\r")
            self.t += float(self.dt)
            FESTIM.solving.solve_it(
                self.F, self.u, J, self.bcs, self.t, self.dt,
                self.parameters["solving"])
            self.u_n.assign(self.u)
            self.post_processing()

    def post_processing(self):
        ret_H, ret_He = 0, 0
        tot_He, tot_H, tot_HeH, tot_V, tot_VH, tot_VHe, tot_VHeH = \
            0, 0, 0, 0, 0, 0, 0

        for cluster in self.clusters:
            if cluster.n_H >= 1:  # lV.mHe.nH
                ret_H += assemble(cluster.n_H*cluster.concentration*dx)
                if cluster.n_He == 0:  # lV.nH
                    if cluster.n_V == 0:  # nH
                        tot_H += cluster.concentration
                    else:  # lV.nH
                        tot_VH += cluster.concentration
                else:  # lV.mHe.nH
                    if cluster.n_V == 0:  # mHe.nH
                        tot_HeH += cluster.concentration
                    else:  # lV.mHe.nH
                        tot_VHeH += cluster.concentration
            elif cluster.n_He >= 1:  # lV.mHe
                ret_He += assemble(cluster.n_He*cluster.concentration*dx)
                if cluster.n_V == 0:  # mHe
                    tot_He += cluster.concentration
                else:  # lV.mHe
                    tot_VHe += cluster.concentration
            elif cluster.n_V >= 1:  # lV
                tot_V += cluster.concentration
        files = [self.file_He_tot, self.file_H_tot, self.file_V_tot,
                 self.file_VH_tot, self.file_VHe_tot, self.file_VHeH_tot]
        fields = [tot_He, tot_H, tot_V, tot_VH, tot_VHe, tot_VHeH]
        names = ["tot_He", "tot_H", "tot_V", "tot_VH", "tot_VHe", "tot_VHeH"]
        self.file_out.write(self.u, self.t)
        for f, field, name in zip(files, fields, names):
            W = self.u.function_space().sub(0).collapse()
            if field == 0:
                field = Constant(0)
            field = project(field, W)
            field.rename(name, "solution")
            f.write(field, self.t)

        self.derived_quantities.append([self.t, self.T, ret_H, ret_He])


if __name__ == "__main__":
    size = 10e-3
    t_final = 1000*3600
    mesh_parameters = {
        "size": size,
        "initial_number_of_cells": 900,
        "refinements": [
            {
                "x": 2e-5,
                "cells": 50
            },
            {
                "x": 2e-6,
                "cells": 600
            },
            {
                "x": 10-9,
                "cells": 600
            },
        ]
    }
    solving_parameters = {
        "adaptive_stepsize": {
            "t_stop": t_final*2,
            "stepsize_stop_max": 10,
            "stepsize_change_ratio": 1.2,
            "dt_min": 1e-5
        },
        "initial_stepsize": 0.5,
        "final_time": t_final,
        "newton_solver": {
            "absolute_tolerance": 1e10,
            "relative_tolerance": 1e-10,
            "maximum_iterations": 15
        },
    }
    flux_He, flux_H = 1e15, 1e15  # /m2/s
    center = 1.5e-9
    width = 1e-9
    distribution = 1/(width*(2*3.14)**0.5) * \
        sp.exp(-0.5*((FESTIM.x-center)/width)**2)
    sources = [
        {
            "cluster": {
                "n_H": 0,
                "n_He": 1,
                "n_V": 0,
            },
            "distribution": distribution,
            "value": flux_He
        },
        {
            "cluster": {
                "n_H": 1,
                "n_He": 0,
                "n_V": 0,
            },
            "distribution": distribution,
            "value": flux_H
        }
    ]
    parameters = {
        "solving": solving_parameters,
        "mesh": mesh_parameters,
        "temperature": 300,
        "sources": sources}

    mySim = Simulation(
        n_He_max=5, n_H_max=5, n_Vac_max=1, parameters=parameters)
    mySim.run()
