import sys


def adaptive_stepsize(nb_it, converged, dt, dt_min, stepsize_change_ratio, t):
    """Adapts the stepsize as function of the number of iterations of the
    solver.

    Arguments:
        nb_it {int} -- number of iterations
        converged {bool} -- True if converged, else False
        dt {fenics.Constant()} -- stepsize
        dt_min {float} -- minimum stepsize
        stepsize_change_ratio {float} -- stepsize change ratio for adaptive
            stepsize
        t {float} -- time
    """

    if converged is False:
        dt.assign(float(dt) / stepsize_change_ratio)
        if float(dt) < dt_min:
            sys.exit("Error: stepsize reached minimal value")

    if nb_it < 5:
        dt.assign(float(dt) * stepsize_change_ratio)
    else:
        dt.assign(float(dt) / stepsize_change_ratio)
    return
