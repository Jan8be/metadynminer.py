import pytest
import metadynminer as mm
from matplotlib import pyplot as plt
import os
import numpy as np

def test_1p():
    expected = np.array([ 385.17567651, 1528.33533688])
    #load hills
    h1 = mm.Hills(name="./data/acealanme1d", periodic=[True])
    #find minima on FES
    minima = mm.Minima(mm.Fes(h1, resolution=256, original=False))
    fep = mm.FEProfile(minima, h1)
    fep = fep.feprofile[:,2:]
    fep = fep.sum(axis=0)
    assert(np.allclose(fep, expected, 1e-3))

def test_2p():
    expected = np.array([ 371.9108889 ,  629.40182382, 1252.17936501, 3040.47957346])
    #load hills
    h2 = mm.Hills(name="./data/acealanme", periodic=[True, True])
    #find minima on FES
    minima = mm.Minima(mm.Fes(h2, resolution=256, original=False))
    fep = mm.FEProfile(minima, h2)
    fep = fep.feprofile[:,2:]
    fep = fep.sum(axis=0)
    assert(np.allclose(fep, expected, 1e-3))

def test_3p():
    expected = np.array([  -41.67918722,  -135.41854014,   243.58920429,  1060.11504292,
        2378.11540026,  6513.85008335,  6598.73851877,  6418.38484675,
        7694.86520632,  7869.748866  ,  8185.64504805,  8250.71152005,
        9633.41481942, 12238.44803989])
    #load hills
    h3 = mm.Hills(name="./data/acealanme3d", periodic=[True, True, True])
    #find minima on FES
    minima = mm.Minima(mm.Fes(h3, resolution=64, original=False))
    fep = mm.FEProfile(minima, h3)
    fep = fep.feprofile[:,2:]
    fep = fep.sum(axis=0)
    assert(np.allclose(fep, expected, 1e-3))