from bubble_growth import main
import numpy as np
import sympy as sp
import FESTIM

# tendril
size = 30e-9
mesh_parameters = {
    "size": size,
    "initial_number_of_cells": 50,
    # "refinements": [
    #     {
    #         "x": 2e-6,
    #         "cells": 600
    #     },
    #     {
    #         "x": 500-9,
    #         "cells": 500
    #     },

    # ]
}

flux = 1e22
folder_name = "faney/tendrils/1000K"

center = 1.5e-9
width = 1e-9
distribution = 1/(width*(2*3.14)**0.5) * \
    sp.exp(-0.5*((FESTIM.x-center)/width)**2)
source = distribution*flux/0.93

main(
    mesh_parameters, dt=1, t_final=5e25/flux, temperature=1000,
    source=source, folder=folder_name)
