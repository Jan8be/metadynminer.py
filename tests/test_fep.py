import pytest
import metadynminer as mm
from matplotlib import pyplot as plt
import os
import numpy as np

def test_1p():
    expected = np.array([ 385.17552302, 1528.3351816 ])
    #load hills
    h1 = mm.Hills(name="./data/acealanme1d", periodic=[True])
    #find minima on FES
    minima = mm.Minima(mm.Fes(h1, resolution=256, original=False))
    fep = mm.FEProfile(minima, h1)
    fep = fep.feprofile[:,2:]
    fep = fep.sum(axis=0)
    assert(np.allclose(fep, expected, 1e-3))

def test_2p():
    expected = np.array([ 371.91080136,  629.40164896, 1289.9170257, 1252.1794877, 3040.4783994 ])
    #load hills
    h2 = mm.Hills(name="./data/acealanme", periodic=[True, True])
    #find minima on FES
    minima = mm.Minima(mm.Fes(h2, resolution=256, original=False))
    fep = mm.FEProfile(minima, h2)
    fep = fep.feprofile[:,2:]
    fep = fep.sum(axis=0)
    assert(np.allclose(fep, expected, 1e-3))

def test_3p():
    expected = np.array([  -41.67860026,   243.58946343,  1237.95317909,  1060.11536049,
        2378.11588445,  6418.3849335 ,  6513.84977833,  7869.74886734,
        7694.86571792,  8185.64534757,  8250.71114721,  9633.41434977,
       12238.44830729])
    #load hills
    h3 = mm.Hills(name="./data/acealanme3d", periodic=[True, True, True])
    #find minima on FES
    minima = mm.Minima(mm.Fes(h3, resolution=64, original=False))
    fep = mm.FEProfile(minima, h3)
    fep = fep.feprofile[:,2:]
    fep = fep.sum(axis=0)
    assert(np.allclose(fep, expected, 1e-3))