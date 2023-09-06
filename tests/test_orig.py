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
    metadynminer_cv1_orig = mm.Fes(h1, resolution=256, original=True).fes.T
    #load plumed FES
    plumed1 = np.loadtxt("./tests/plumed_acealanme1d.dat")
    plumed1 = np.reshape(plumed1[:,1], (256))
    plumed1 = plumed1 - np.min(plumed1)
    # compare
    abs_diff = np.max(np.abs(plumed1-metadynminer_cv1_orig))
    assert(abs_diff < 1e-6)
    
def test_2p():
    expected = 1810417.2978640283
    #load hills
    h2 = mm.Hills(name="./data/acealanme", periodic=[True,True])
    #prepare FES
    metadynminer_cv2_orig = mm.Fes(h2, resolution=256, original=True).fes.T
    #load plumed FES
    plumed2 = np.loadtxt("./tests/plumed_acealanme.dat")
    plumed2 = np.reshape(plumed2[:,2], (256,256))
    plumed2 = plumed2 - np.min(plumed2)
    # compare
    abs_diff = np.max(np.abs(plumed2-metadynminer_cv2_orig))
    assert(abs_diff < 1e-6)
    
def test_3p():
    expected = 15149816.57067392
    #load hills
    h3 = mm.Hills(name="./data/acealanme3d", periodic=[True,True,True])
    #prepare FES
    metadynminer_cv3_orig = mm.Fes(h3, resolution=64, original=True).fes.T
    #load plumed FES
    plumed3 = np.loadtxt("./tests/plumed_acealanme3d_64.dat")
    plumed3 = np.reshape(plumed3[:,3], (64,64,64))
    plumed3 = plumed3 - np.min(plumed3)
    # compare
    abs_diff = np.max(np.abs(plumed3-metadynminer_cv3_orig))
    assert(abs_diff < 1e-6)
    
def test_1np():
    expected = 11530.43835598843
    #load hills
    h1 = mm.Hills(name="./data/acealanme1d", periodic=[False])
    #prepare FES
    metadynminer_cv1_fast = mm.Fes(h1, resolution=256, original=True).fes.T
    #sum FES
    mm_sum = np.sum(metadynminer_cv1_fast)
    #compare
    abs_diff = np.abs(expected - mm_sum)
    assert(abs_diff < 1e-4)

def test_2np():
    expected = 3541846.031986131
    #load hills
    h2 = mm.Hills(name="./data/acealanme", periodic=[False,False])
    #prepare FES
    metadynminer_cv2_fast = mm.Fes(h2, resolution=256, original=True).fes.T
    #sum FES
    mm_sum = np.sum(metadynminer_cv2_fast)
    #compare
    abs_diff = np.abs(expected - mm_sum)
    assert(abs_diff < 1e-4)
    
def test_3np():
    expected = 14173533.11781302
    #load hills
    h3 = mm.Hills(name="./data/acealanme3d", periodic=[False,False,False])
    #prepare FES
    metadynminer_cv3_fast = mm.Fes(h3, resolution=64, original=True).fes.T
    #sum FES
    mm_sum = np.sum(metadynminer_cv3_fast)
    #compare
    abs_diff = np.abs(expected - mm_sum)
    assert(abs_diff < 1e-4)
    
if __name__ == '__main__':
    pytest.main([__file__])
