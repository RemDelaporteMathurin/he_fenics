from bubble_growth import main
import numpy as np
import sympy as sp
import FESTIM

nb_points = 20
min_T = 100
max_T = 1200
temperatures = np.random.uniform(low=min_T, high=max_T, size=(nb_points,))

min_flux = 17+np.log10(1)
max_flux = 21+np.log10(5)
fluxes = np.random.uniform(low=min_flux, high=max_flux, size=(nb_points,))
points = np.c_[temperatures, fluxes]
size = 100e-6
mesh_parameters = {
    "size": size,
    "initial_number_of_cells": 400,
    "refinements": [
        {
            "x": 2e-6,
            "cells": 600
        },
        {
            "x": 50-9,
            "cells": 50
        },
        {
            "x": 8e-9,
            "cells": 25
        },
    ]
}
if __name__ == "__main__":
    for point in points:
        flux = 10**point[1]
        folder_name = \
            "parametric_study/T={:.3e}_flux={:.2e}".format(point[0], flux)

        center = 1.5e-9
        width = 1e-9
        distribution = 1/(width*(2*3.14)**0.5) * \
            sp.exp(-0.5*((FESTIM.x-center)/width)**2)
        source = distribution*flux/0.93

        main(
            mesh_parameters, dt=0.000001, t_final=3600, temperature=point[0],
            source=source, folder=folder_name)
