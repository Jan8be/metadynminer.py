import pytest
import metadynminer as mm
from matplotlib import pyplot as plt
import os
import numpy as np

def test_it():
  #load hills
  h1 = mm.Hills(name="acealanme1d", periodic=[True])
  #prepare FES
  metadynminer_cv1_fast = mm.Fes(h1, resolution=256, original=False).fes.T
  #load plumed FES
  plumed1 = np.loadtxt("plumed_acealanme1d.dat")
  plumed1 = np.reshape(plumed1[:,1], (256))
  plumed1 = plumed1 - np.min(plumed1)
  #compare
  mean_error_limit_fast = 1
  max_error_limit_fast = 4
  mean_error = np.mean(metadynminer_cv1_fast-plumed1)
  max_error = np.max(np.abs(metadynminer_cv1_fast-plumed1))
  assert(max_error < max_error_limit_fast & mean_error < mean_error_limit_fast)

if __name__ == '__main__':
  pytest.main([__file__])
