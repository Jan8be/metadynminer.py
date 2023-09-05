import pytest
import metadynminer as mm
from matplotlib import pyplot as plt
import os
import numpy as np

def test_1p():
    expected = 4345.243401728701
    #load hills
    h1 = mm.Hills(name="./data/acealanme1d", periodic=[True])
    #prepare FES
    metadynminer_cv1_fast = mm.Fes(h1, resolution=256, original=False).fes.T
    #sum FES
    mm_sum = np.sum(metadynminer_cv1_fast)
    #compare
    abs_diff = np.abs(expected - mm_sum)
    assert(abs_diff < 1e-4)
    
def test_1np():
    expected = 11586.547055427745
    #load hills
    h1 = mm.Hills(name="./data/acealanme1d", periodic=[False])
    #prepare FES
    metadynminer_cv1_fast = mm.Fes(h1, resolution=256, original=False).fes.T
    #sum FES
    mm_sum = np.sum(metadynminer_cv1_fast)
    #compare
    abs_diff = np.abs(expected - mm_sum)
    assert(abs_diff < 1e-4)
    
def test_2p():
    expected = 1810417.2978640283
    #load hills
    h2 = mm.Hills(name="./data/acealanme", periodic=[True,True])
    #prepare FES
    metadynminer_cv2_fast = mm.Fes(h2, resolution=256, original=False).fes.T
    #sum FES
    mm_sum = np.sum(metadynminer_cv2_fast)
    #compare
    abs_diff = np.abs(expected - mm_sum)
    assert(abs_diff < 1e-4)

def test_2np():
    expected = 3569701.989450706
    #load hills
    h2 = mm.Hills(name="./data/acealanme", periodic=[False,False])
    #prepare FES
    metadynminer_cv2_fast = mm.Fes(h2, resolution=256, original=False).fes.T
    #sum FES
    mm_sum = np.sum(metadynminer_cv2_fast)
    #compare
    abs_diff = np.abs(expected - mm_sum)
    assert(abs_diff < 1e-4)
    
def test_3p():
    expected = 15149816.57067392
    #load hills
    h3 = mm.Hills(name="./data/acealanme3d", periodic=[True,True,True])
    #prepare FES
    metadynminer_cv3_fast = mm.Fes(h3, resolution=64, original=False).fes.T
    #sum FES
    mm_sum = np.sum(metadynminer_cv3_fast)
    #compare
    abs_diff = np.abs(expected - mm_sum)
    assert(abs_diff < 1e-4)

#def test_3np():
#    expected = 15149816.57067392
#    #load hills
#    h3 = mm.Hills(name="./data/acealanme3d", periodic=[False,False,False])
#    #prepare FES
#    metadynminer_cv3_fast = mm.Fes(h3, resolution=64, original=False).fes.T
#    #sum FES
#    mm_sum = np.sum(metadynminer_cv3_fast)
#    #compare
#    abs_diff = np.abs(expected - mm_sum)
#    assert(abs_diff < 1e-4)

if __name__ == '__main__':
    pytest.main([__file__])
