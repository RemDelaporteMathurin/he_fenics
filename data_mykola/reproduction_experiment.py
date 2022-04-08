import os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

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
            "cells": 50
        },
    ]
}

flux = 2.3e22
folder_name = "simulation"

center = 1.5e-9
width = 0.8e-9
distribution = 1/(width*(2*3.14)**0.5) * \
    sp.exp(-0.5*((FESTIM.x-center)/width)**2)
source = distribution*flux*0.2

main(
    mesh_parameters, dt=1e-5, t_final=3e23/flux, temperature=1053,
    source=source, folder=folder_name)
