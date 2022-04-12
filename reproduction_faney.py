from bubble_growth import main
import numpy as np
import sympy as sp
import FESTIM
from FESTIM import x, t

# tendril
size = 30e-9
mesh_parameters = {
    "size": size,
    "initial_number_of_cells": 300,
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
#folder_name = "faney/tendrilsJ/1000K"
folder_name = ""

# TS#1
#center = 1.5e-9
#width = 1e-9
#distribution = 1/(width*(2*3.14)**0.5) * sp.exp(-0.5*((FESTIM.x-center)/width)**2)
#print('Integration of distribution:', float(sp.integrate(distribution, (x, 0.1, 30e-9))))
#distribution = distribution/sp.integrate(distribution, (x, 0, 30e-9))

# TS Xolotl PSI
#params= [-0.00200307652347,0.045491357931,0.257229692641,-0.431991663267,0.374291001257,-0.205933999597,0.0769701218727,-0.020277631463,0.00384532323736,-0.000530022344926,5.30860448111e-05,-3.81962426235e-06,1.92207763422e-07,-6.41748946396e-09,1.27665777741e-10,-1.1449383864e-12]
#totalDepths = 14.6512438794
#distribution = (x <= totalDepths*1e-9) * 1e9 * (params[0] + params[1] * x*1e9 + params[2] * pow(x*1e9, 2.0) + params[3] * pow(x*1e9, 3.0) + params[4] * pow(x*1e9, 4.0) + params[5] * pow(x*1e9, 5.0) + params[6] * pow(x*1e9, 6.0) + params[7] * pow(x*1e9, 7.0) + params[8] * pow(x*1e9, 8.0) + params[9] * pow(x*1e9, 9.0) + params[10] * pow(x*1e9, 10.0) + params[11] * pow(x*1e9, 11.0) + params[12] * pow(x*1e9, 12.0) + params[13] * pow(x*1e9, 13.0) + params[14] * pow(x*1e9, 14.0) + params[15] * pow(x*1e9, 15.0))

# TS Xolotl W100
distribution = (x <= 10e-9) *  1e9/16.39 *  (7.00876507 + 0.6052078 * x*1e9 - 3.01711048 *(x*1e9)**2 + 1.36595786 * (x*1e9)**3 - 0.295595 * (x*1e9)**4 + 0.03597462 * (x*1e9)**5 - 0.0025142 * (x*1e9)**6 + 0.0000942235 * (x*1e9)**7 - 0.0000014679 * (x*1e9)**8)


source = distribution*flux

main(
    mesh_parameters, dt=1000.0, t_final=1e12, temperature=1500,
    source=source, folder=folder_name)
