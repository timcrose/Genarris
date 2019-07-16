# This file is part of aimsutils
# (C) 2015 Christoph Schober

import sys
from collections import OrderedDict

import numpy as np
from numpy import sqrt as sqrt2


def construct_np_C(l):
    """
    Construct the C matrix according to Blanco et al (eq. 19).
    This version uses numpy instead sympy and is therefore considerably
    faster, although limited to numpy precision (which usually is not a
    problem).

    Parameters
    ----------
    l : int
        Angular momentum l.

    Returns
    -------
    C : np.array
        The C matrix to transform complex to real YLMs.
    """
    ms = np.arange(-l, l+1, 1)

    C = np.zeros([len(ms), len(ms)], dtype=complex)

    for m in ms:
        m = int(m)
        for m2 in ms:
            m2 = int(m2)
            if abs(m) != abs(m2):
                element = 0
            elif m == 0 and m2 == 0:
                element = sqrt2(2)
            elif abs(m) > 0:
                if np.sign(m) == np.sign(m2) == -1:
                    # upper left, m<0 + m2<0
                    element = 1j
                elif np.sign(m) == 1 and np.sign(m2) == -1:
                    # lower left, m>0 + m2<0
                    element = 1
                elif np.sign(m) == 1 == np.sign(m2) == 1:
                    # lower right, m>0 + m2>0
                    element = (-1)**(l-(l-m))
                elif np.sign(m) == -1 and np.sign(m2) == 1:
                    # upper right, m<0 + m2 > 0
                    element = -1j*(-1)**(l-(l-m))

            C[m+l, m2+l] = element
    C = (1/sqrt2(2))*C

    return C
