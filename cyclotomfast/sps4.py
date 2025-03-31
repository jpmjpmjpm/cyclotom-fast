from copy import deepcopy

from sympy.functions.combinatorial.numbers import totient
from sympy.ntheory import factorint
# from utils.debug import dump_args
import numpy as np


class SqFreeFactors:
    """
    A class used to represent square free factors of an integer

    value: integer value
    totient: Euler's totient of the integer
    factors: the prime factors of the integer
    """
    value = None
    totient = None
    factors = None

    def __init__(self, value: int) -> None:
        self.value = value
        self.totient = totient(value)
        dict_factors = factorint(value)
        self.factors = [*dict_factors]
        self.factors.sort()

    def __repr__(self):
        return f"{self.factors}"

    def has_prime(self) -> bool:
        return self.value > 1

    def remove_last(self):
        prime = self.factors[-1]
        self.value //= prime
        self.totient //= prime - 1
        del self.factors[-1]
        return prime


# =======================================================================================
# NUMPY versions of polynomial operations
# =======================================================================================
def np_polymul(poly, n):
    """Multiply by 1 -x^n """
    poly0 = [0] * (n + 1)
    poly0[0], poly0[-1] = -1, 1
    poly = np.polymul(poly, np.poly1d(poly0))
    return poly


def np_polydiv(poly, n):
    """Divide by 1 -x^n """
    poly0 = [0] * (n + 1)
    poly0[0], poly0[-1] = -1, 1
    poly, _ = np.polydiv(poly, np.poly1d(poly0))
    return poly


def np_first(alpha: bool):
    return np.poly1d([1 if alpha else -1])


# =======================================================================================
# SERIES versions of polynomial operations
# =======================================================================================
class Polyc:
    """
    degree: degree of the polynomial
    coeffs: coefficients of the polynomial
    """
    deg: None
    coeffs: []

    def __init__(self, deg=None, coeffs=None):
        self.deg = deg if deg else len(coeffs) - 1
        self.coeffs = coeffs.copy() if coeffs else [0] * (deg + 1)

    @staticmethod
    def change_dim(polyc, dim):
        new_polyc = Polyc(dim)
        new_polyc.coeffs = [polyc.coeffs[i] for i in range(max(polyc.deg, dim))]
        return new_polyc


def sr_polymul(poly, n):
    """ Multiply by 1 -x^n
        TODO Limit number of copies """
    coeffs = poly.coeffs.copy() + [0] * n
    for i in range(poly.deg, -1, -1):
        coeffs[i + n] -= coeffs[i]
    return Polyc(coeffs=coeffs)


def sr_polydiv(poly, n):
    """ Divide by 1 - x^n
        TODO Look at why only 1 + x^n is computed in the paper """
    deg = poly.deg
    degd = deg - n
    coeffsd = poly.coeffs[:-n].copy()
    for power in range(n, degd + 1, n):
        for i in range(0, degd - power + 1):
            coeffsd[i + power] += poly.coeffs[i]
    return Polyc(coeffs=coeffsd)


def sr_first(alpha: bool):
    return Polyc(coeffs=[1 if alpha else -1])


def cyclotomic(m: SqFreeFactors, alpha: bool, polyfirst, polymul, polydiv):
    """Return a cyclotomic or inverse cyclotomic polynomial.

    Parameters
    ==========
    m : SqFreeFactors
        Square free odd integer
    alpha: bool
        True to compute the cyclotomic polynomial, otherwise compute the inverse cyclotomic polynomial
    polyfirst:
        Method to compute the initial polynomial: +1 for cyclotomic or -1 for inverse cyclotomic
    polymul:
        Method to compute multiplication by 1 - z^n
    polydiv:
        Method to compute division by 1 - z^n
    """
    return sps4(m, 1, True, polyfirst(alpha), polymul, polydiv)


# @dump_args
def sps4(m: SqFreeFactors, e: int, alpha: bool, poly, polymul, polydiv):
    """
    Multiply poly polynomial by Phi_m(z^e) if alpha is True otherwise by Psi_m(z^e)
    :param m: SqFreeFactors
        Square free odd integer
    :param e: int
        The power of z in the polynomial used for the product
    :param alpha: bool
        Multiply by Phi_m(z^e) if alpha is True otherwise by Psi_m(z^e)
    :param poly: the polynomial to be multiplied
    :param polymul: method to compute multiplication by 1 - z^n
    :param polydiv: method to compute division by 1- z^n
    :return: the product of the polynomials
    """
    ms, es = deepcopy(m), e
    while ms.has_prime():
        p = ms.remove_last()
        poly = sps4(ms, es, not alpha, poly, polymul, polydiv)
        es *= p

    if alpha:
        # Multiply by 1 - z^{m·e}
        poly = polymul(poly, m.value * e)

        # For each prime factor p of m
        for p in m.factors.copy():
            # Divide by 1 - z^{m·e/p}
            poly = polydiv(poly, m.value * e // p)

    return poly
