from copy import deepcopy
from typing import List
from sympy.functions.combinatorial.numbers import totient
from sympy.ntheory import factorint
from utils.tests import dump_args
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


@dump_args
def sps4(m: SqFreeFactors, e: int, alpha: bool, poly):
    """
    Multiply a polynomial by \Phi_m(z^e) if alpha is True otherwise by \Psi_m(z^e)
    :param m: SqFreeFactors
        Square free odd integer
    :param e: int
        The power of z in the polynomial used for the product
    :param alpha: bool
        Multiply by \Phi_m(z^e) if alpha is True otherwise by \Psi_m(z^e)
    :param poly: the polynomial to be multiplied
    :return: the product of the polynomials
    """
    ms, es = deepcopy(m), e
    while ms.has_prime():
        p = ms.remove_last()
        poly = sps4(ms, es, not alpha, poly)
        es *= p

    if alpha:
        # Multiply by 1 - z^{m·e}
        poly0 = [0] * (m.value * e + 1)
        poly0[0], poly0[-1] = -1, 1
        poly = np.polymul(poly, np.poly1d(poly0))

        # For each prime factor p of m
        for p in m.factors.copy():
            # Divide by 1 - z^{m·e/p}
            poly0 = [0] * (m.value * e // p + 1)
            poly0[0], poly0[-1] = -1, 1
            poly, _ = np.polydiv(poly, np.poly1d(poly0))

    return poly


def cyclotomic(m: SqFreeFactors, alpha: bool):
    """Return a cyclotomic or inverse cyclotomic polynomial.

    Parameters
    ==========
    m : SqFreeFactors
        Square free odd integer
    alpha: bool
        True to compute the cyclotomic polynomial, otherwise compute the inverse cyclotomic polynomial
    """
    constant_term = 1 if alpha else -1
    poly = np.poly1d([constant_term])
    return sps4(m, 1, True, poly)
