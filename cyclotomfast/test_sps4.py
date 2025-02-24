from sps4 import cyclotomic, SqFreeFactors, np_first, np_polymul, np_polydiv
import numpy as np


def test_cyclotomic_5():
    sqf = SqFreeFactors(5)
    assert cyclotomic(sqf, True, np_first, np_polymul, np_polydiv) == np.poly1d([1., 1., 1., 1., 1.])


def test_cyclotomic_15():
    sqf = SqFreeFactors(15)
    assert cyclotomic(sqf, True, np_first, np_polymul, np_polydiv) == np.poly1d(
        [1., -1., 0., 1., -1., 1.0, 0., -1., 1.0])


def test_cyclotomic_105():
    sqf_105 = SqFreeFactors(105)
    phi_105 = cyclotomic(sqf_105, True, np_first, np_polymul, np_polydiv)
    assert phi_105[7] == -2.0
