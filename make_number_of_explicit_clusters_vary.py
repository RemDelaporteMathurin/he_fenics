from bubble_growth import main
import numpy as np
import sympy as sp
import FESTIM

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

    flux = 1e20
    center = 1.5e-9
    width = 1e-9
    distribution = 1/(width*(2*3.14)**0.5) * \
        sp.exp(-0.5*((FESTIM.x-center)/width)**2)
    source = distribution*flux/0.93
    for N in [20]:
        main(
            mesh_parameters, dt=0.000001, t_final=50, temperature=1000,
            source=source, nb_clusters=N, folder="make_N_vary/N={}".format(N))
