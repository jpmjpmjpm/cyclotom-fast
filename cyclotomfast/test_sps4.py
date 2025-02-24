from sps4 import cyclotomic, SqFreeFactors
import numpy as np


def test_cyclotomic_5():
    sqf = SqFreeFactors(5)
    assert cyclotomic(sqf, True) == np.poly1d([1., 1., 1., 1., 1.])


def test_cyclotomic_15():
    sqf = SqFreeFactors(15)
    assert cyclotomic(sqf, True) == np.poly1d([1., -1., 0., 1., -1., 1.0, 0., -1., 1.0])


def test_cyclotomic_105():
    sqf_105 = SqFreeFactors(105)
    phi_105 = cyclotomic(sqf_105, True)
    assert phi_105[7] == -2.0
