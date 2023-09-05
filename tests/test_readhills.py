import pytest
import metadynminer as mm
from matplotlib import pyplot as plt
import os
import numpy as np

def test_1p():
    expected = 1025.0783977216151
    #load hills
    h1 = mm.Hills(name="./data/acealanme1d", periodic=[True])
    #sum hills heights in HILLS
    hills_sum = np.sum(h1.heights)
    assert(np.allclose(hills_sum, expected))
    
def test_2p():
    expected = 4861.35454589019
    #load hills
    h2 = mm.Hills(name="./data/acealanme", periodic=[True,True])
    #sum hills heights in HILLS
    hills_sum = np.sum(h2.heights)
    assert(np.allclose(hills_sum, expected))
    
def test_3p():
    expected = 11610.060496549602
    #load hills
    h3 = mm.Hills(name="./data/acealanme3d", periodic=[True,True,True])
    #sum hills heights in HILLS
    hills_sum = np.sum(h3.heights)
    assert(np.allclose(hills_sum, expected))

if __name__ == '__main__':
    pytest.main([__file__])
