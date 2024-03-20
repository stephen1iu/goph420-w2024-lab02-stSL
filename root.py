import numpy as np

def root_newton_raphson(x0, f, dfdx):
    """
    Inputs
    ------
    x0: initial guess
    f: function
    dfdx: first order derivative of f

    Returns
    -------
    float: estimate of root
    int: number of iterations for convergence
    numpy.ndarray: approximate relative errors
    """

    x=x0
    i=0
    eps_a=[1]
    tol=1e-8

    while np.abs(eps_a[i])>tol:
        x_old=x
        x=(x-(f(x)/dfdx(x)))
        eps_a.append((x-x_old)/(x))
        i+=1
        if i>150:
            raise RuntimeWarning ("no root found")
    return(
        x,
        i,
        eps_a[1:] #doesnt return intializer value of eps_a
    )