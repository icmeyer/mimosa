import numpy as np
import cmath
import re

def checktol(x,y,tol):
    """Check absolute difference between two values and compare to 
    a defined tolerance. Return boolean. """
    err = abs(abs(x-y)/x)
    if err>tol:
        return False
    elif err<tol:
        return True

def collect_phi(regions, adjoint=False):
    """
    Arrange the phi of all regions into a single vector
    """
    ngroup = regions[0].phi.shape[0]
    phi = np.zeros(len(regions)*ngroup)
    counter = 0
    for region in regions:
        if adjoint:
            phi[ngroup*counter:ngroup*(counter+1)] = region.a_phi
        else:
            phi[ngroup*counter:ngroup*(counter+1)] = region.phi
    return phi
        

def check_phi_convergence(regions, old_phi, adjoint=False):
    """
    Check for convergence across phi using old_vector
    """
    new_phi = collect_phi(regions, adjoint=adjoint)
    err = np.linalg.norm(new_phi - old_phi)/np.linalg.norm(old_phi)
    # print(err)
    if err < 1e-6:
        converge_flag = True
    else:
        converge_flag = False

    return converge_flag, new_phi

def normalize(x):
    return x/np.linalg.norm(x)

def get_trailing_numbers(s, zero=False):
    m = re.search(r'\d+$', s)
    if m:
        return int(m.group())
    elif zero:
        return 0
    else:
        return None

def real_reduce(num):
    """ Reduce complex numbers to reals if imaginary part is 0 """
    return num.real if (num.imag == 0).all() else num

def quad_sum(array):
    print(array)
    n = len(array)
    total = 0
    for i in range(n):
        for j in range(n):
            total += array[i,j]**2
    return total
