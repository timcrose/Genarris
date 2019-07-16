# This file is part of aimsutils.
# (C) 2015 Christoph Schober

import numpy as np
from collections import OrderedDict
import math
import sys

from sympy import Matrix, Symbol, cos, sin, pi as PI, sqrt as SQRT
from aimsutils.rotate import wigner


def ylm(theta, phi, lmin, lmax):
    """
    Symbolic Python version of FHI aims YLM algorithm (public domain). The
    function uses sympy for the symbolic equations. This python version is
    identical with the AIMS YLMs in ``increment_ylm.f90``.

    Parameters
    ----------
    theta : float
        Angle theta in radians.
    phi : float
        Angle phi in radians.
    lmin : int
        Lowest angular momentum l to calculate.
    lmax : int
        Maximum angular momentum l to calculate.

    Returns
    -------
    sY : OrderedDict
        Dictionary with all YLMs for given parameters (theta, phi, lmin, lmax).
    """

    COSTH = cos(theta)
    SINTH = sin(theta)

    COSPH = cos(phi)
    SINPH = sin(phi)

    Y = [0 for x in range((lmax+1)*(lmax+1)+1)]
    if lmin == 0:
        YLLR = 1./SQRT(4.*PI)
        YLLI = 0
        Y[1] = YLLR

    if lmin <= 1 and lmax >= 1:
        Y[3] = SQRT(3.)*YLLR*COSTH

        TEMP1 = -SQRT(3.)*YLLR*SINTH
        Y[4] = TEMP1*COSPH
        Y[2] = -TEMP1*SINPH

    for L in np.arange(max(2, lmin), lmax+1):
        INDEX = L*L+1
        INDEX2 = INDEX + 2*L
        MSIGN = 1 - 2*np.mod(L, 2)

        YL1L1R = Y[INDEX-1]
        YL1L1I = -MSIGN * Y[INDEX-2*L+1]
        TEMP1 = -SQRT((2.*L+1)/(2.*L))*SINTH
        YLLR = TEMP1*(COSPH*YL1L1R - SINPH*YL1L1I)
        YLLI = TEMP1*(COSPH*YL1L1I + SINPH*YL1L1R)
        Y[INDEX2] = YLLR
        Y[INDEX] = MSIGN * YLLI
        INDEX2 = INDEX2 - 1
        INDEX = INDEX + 1

        TEMP2 = SQRT((2.*L+1))*COSTH
        YLL1R = TEMP2*YL1L1R
        YLL1I = TEMP2*YL1L1I
        Y[INDEX2] = YLL1R
        Y[INDEX] = -MSIGN * YLL1I
        INDEX2 = INDEX2 - 1
        INDEX = INDEX + 1

        I4L2 = INDEX - 4*L + 2
        I2L = INDEX - 2*L
        I24L2 = INDEX2 - 4*L + 2
        I22L = INDEX2 - 2*L
        D4LL1C = COSTH*SQRT((4.*L*L-1))
        D2L13 = -SQRT((2.*L+1)/(2.*L-3))

        for M in range(L-2, -1, -1):
            TEMP1 = 1./SQRT(((L+M)*(L-M)))
            TEMP2 = D4LL1C*TEMP1
            TEMP3 = D2L13*SQRT(((L+M-1)*(L-M-1)))*TEMP1
            YLMR = TEMP2*Y[I22L] + TEMP3*Y[I24L2]
            YLMI = TEMP2*Y[I2L] + TEMP3*Y[I4L2]
            Y[INDEX2] = YLMR
            Y[INDEX] = YLMI

            INDEX2 = INDEX2 - 1
            INDEX = INDEX + 1
            I24L2 = I24L2 - 1
            I22L = I22L - 1
            I4L2 = I4L2 + 1
            I2L = I2L + 1

    YY = [yy for yy in Y[1:]]

    sY = OrderedDict()
    i = 0
    for l in range(lmax+1):
        m = 2*l+1
        sY[l] = YY[i:i+m]
        i = i+m
    return sY


def flatten_array(in_array):
    """
    Take an array with subarrays for each group of basis functions and return
    a flattened version.

    Parameters
    ----------
    in_array : np.array
        The input array.
    n_basis : int
        Number of basis functions
    n_states : int
        Number of states

    Returns
    -------
    array : np.array
        The flattened version of the input array.
    """
    n_basis = sum([x.shape[0] for x in in_array])
    shape = in_array[0].shape
    n_states = shape[1]
    n_spin = shape[2]

    array = list()
    for i in in_array:
        if len(i.shape) != 1:
            array.extend(i.flatten())
        else:
            array.extend(i)

    array = np.array(array)
    array = np.array(array.reshape(n_basis, n_states, n_spin))

    return array


def group_array(in_array, lvec):
    """
    Take a flattened array and impose the grouping of the basis functions.

    Parameters
    ----------
    in_array : np.array
        The input array.
    lvec : list
        A list with all angular momenta.

    Returns
    -------
    array : np.array
        The input array with batches per basis function.
    """

    array = []

    line = 0
    while line < len(lvec):
        l = lvec[line]
        m = 2*l+1
        array.append(np.array(in_array[line:line+m]))
        line += m

    return array


def get_l_vector(basis):
    """
    Construct vector with all l in AIMS order from basis meta data as obtained
    via aimsutils.parser.parse_basis.
    """
    lvec = list()
    for key, value in basis.items():
        lvec.append(value[2])
    return lvec


def aimsify_C(lmax):
    """
    Transform textbook C matrix to aims C matrix.
    """
    C_aims = OrderedDict()
    for l in np.arange(0, lmax+1):
        C_tb = wigner.construct_C(l)
        T = get_num_T(l)[l]
        C_aims[l] = Matrix(np.diag(T))*C_tb

    return C_aims


def get_num_T(lmax):
    """
    Get T matrix by numeric comparison for sympy YLMs with aims YLMs.
    The routine checks for special cases when the value is 0 and tries
    different angles until a non-zero value is found.

    Parameters
    ----------
    lmax : int
        The maximum angular momentum l.

    Returns
    -------
    T : OrderedDict
        A dictionary with the diagonals of the transformation matrix T for
        each l.
    """
    theta = Symbol("theta")
    phi = Symbol("phi")

    aims = ylm(theta, phi, 0, lmax)
    sym = wigner.znm(theta, phi, lmax)

    T = OrderedDict()
    increment = 3
    for l, value in aims.items():
        T[l] = list()
        for m, exp in enumerate(value):
            aims_value = 0
            phi_s = 140
            theta_s = 140
            i = 0
            while abs(aims_value) < 1E-5:
                if i > 0:
                    print("i-th round of while loop")
                sublist = [(phi, math.radians(phi_s+i*increment)),
                           (theta, math.radians(theta_s+i*increment))]
                aims_value = exp.subs(sublist).evalf()
                try:
                    sym_value = (sym[l][m].subs(sublist).evalf())\
                        .as_real_imag()[0]
                except AttributeError:
                    sym_value = (sym[l][m].subs(sublist).evalf())
                i += 1
                test = abs(aims_value) - abs(sym_value)
                if test > 1E-5:
                    sys.exit("Found different value for l:{0} m:{1} ({2})"
                             .format(l, m, test))
                i += 1
            sign = np.sign(aims_value / sym_value)
            T[l].append(sign)

    return T
