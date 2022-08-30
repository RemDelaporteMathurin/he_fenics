from bubble_growth import main
import sympy as sp

x = sp.Symbol("x[0]")
size = 100e-6

mesh_parameters = {
    "size": size,
    "initial_number_of_cells": 900,
    "refinements": [
        {"x": 2e-6, "cells": 600},
        {"x": 50 - 9, "cells": 50},
    ],
}
center = 1.5e-9
width = 1e-9
distribution = (
    1 / (width * (2 * 3.14) ** 0.5) * sp.exp(-0.5 * ((x - center) / width) ** 2)
)
flux = 1e22  # flux in He m-2 s-1
source = distribution * flux / 0.93

dt = 0.0001
t_final = 50
temperature = 1000

# run with bursting
main(
    mesh_parameters,
    dt=dt,
    t_final=t_final,
    temperature=temperature,
    source=source,
    folder="vacancies/bursting",
    k_burst_0=2e3,
)


# run without bursting
main(
    mesh_parameters,
    dt=dt,
    t_final=t_final,
    temperature=temperature,
    source=source,
    folder="vacancies/no_bursting",
    k_burst_0=0,
)