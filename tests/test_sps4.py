from cyclotomfast.sps4 import  (cyclotomic, SqFreeFactors, np_first, np_polymul, np_polydiv,
                  Polyc, sr_first, sr_polymul, sr_polydiv)
import numpy as np

def test_sr_polymul():
    assert sr_polymul(Polyc(coeffs=[1, 2]), 1).coeffs == [1, 1, -2]
    assert sr_polymul(Polyc(coeffs=[1, 2]), 5).coeffs == [1, 2, 0, 0, 0, -1, -2]
    assert sr_polymul(Polyc(coeffs=[1, 2]), 2).coeffs == [1, 2, -1, -2]

def test_sr_polydiv():
    assert sr_polydiv(Polyc(coeffs=[1, 1, -2]), 1).coeffs == [1, 2]
    assert sr_polydiv(Polyc(coeffs=[1, 2, 0, 0, 0, -1, -2]), 5).coeffs == [1, 2]
    assert sr_polydiv(Polyc(coeffs=[1, 2, -1, -2]), 2).coeffs == [1, 2]


def test_np_cyclotomic_5():
    sqf = SqFreeFactors(5)
    assert cyclotomic(sqf, True, np_first, np_polymul, np_polydiv) == np.poly1d([1., 1., 1., 1., 1.])


def test_np_cyclotomic_15():
    sqf = SqFreeFactors(15)
    assert cyclotomic(sqf, True, np_first, np_polymul, np_polydiv) == np.poly1d(
        [1., -1., 0., 1., -1., 1.0, 0., -1., 1.0])


def test_sr_cyclotomic_15():
    sqf = SqFreeFactors(15)
    assert (cyclotomic(sqf, True, sr_first, sr_polymul, sr_polydiv)).coeffs == [1, -1, 0, 1, -1, 1, 0, -1, 1]


def test_np_cyclotomic_105():
    sqf = SqFreeFactors(105)
    phi = cyclotomic(sqf, True, np_first, np_polymul, np_polydiv)
    assert phi[7] == -2.0


def test_sr_cyclotomic_255255():
    sqf = SqFreeFactors(255255)
    phi = cyclotomic(sqf, True, sr_first, sr_polymul, sr_polydiv)
    assert phi.height() == 532


def test_sr_cyclotomic_1181895():
    sqf = SqFreeFactors(1181895)
    phi = cyclotomic(sqf, True, sr_first, sr_polymul, sr_polydiv)
    assert phi.height() == 14102773
