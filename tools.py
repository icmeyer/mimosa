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
