import string

import numpy as np

from .base_constituent import BaseConstituent
from ..nodal_corrections import u_zero, f_unity


class Constituent(BaseConstituent):
    """
    An elementary constituent is fundamental, having its own individual u
    and f functions.
    """
    def __init__(self, name, xdo, u=u_zero, f=f_unity):
        super(Constituent, self).__init__(name)
        self.coefficients = np.array(xdo_to_coefficients(xdo))
        self.u_function = u
        self.f_function = f

    def speed(self, a):
        return np.dot(self.coefficients, self.astro_speeds(a))

    def V(self, astro):
        return np.dot(self.coefficients, self.astro_values(astro))

    def u(self, astro):
        return self.u_function(astro)

    def f(self, astro):
        return self.f_function(astro)

    def xdo(self):
        return coefficients_to_xdo(self.coefficients)


_XDO_TO_INT = {
    'A': 1,
    'B': 2,
    'C': 3,
    'D': 4,
    'E': 5,
    'F': 6,
    'G': 7,
    'H': 8,
    'I': 9,
    'J': 10,
    'K': 11,
    'L': 12,
    'M': 13,
    'N': 14,
    'O': 15,
    'P': 16,
    'Q': 17,
    'R': -8,
    'S': -7,
    'T': -6,
    'U': -5,
    'V': -4,
    'W': -3,
    'X': -2,
    'Y': -1,
    'Z': 0
}

_INT_TO_XDO = {v: k for k, v in _XDO_TO_INT.items()}


def xdo_to_coefficients(xdos):
    return [_XDO_TO_INT[letter.upper()]
            for letter in xdos
            if letter in string.ascii_letters]


def coefficients_to_xdo(coefficients):
    return ''.join([_INT_TO_XDO[c] for c in coefficients])
