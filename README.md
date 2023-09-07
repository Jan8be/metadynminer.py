# metadynminer.py

[![Build](https://github.com/Jan8be/metadynminer.py/actions/workflows/ci.yml/badge.svg)](https://github.com/Jan8be/metadynminer.py/actions/workflows/ci.yml)

Metadynminer is a package designed to help you analyse output HILLS files from PLUMED metadynamics simulations. 

It is based on Metadynminer package for R programming language, but it is not just a port from R to Python, as it is updated and improved in many aspects. It supports HILLS files with one, two or three collective variables. 

All these functions can be easily customized with many parameters. You can learn more about that later in the documentation. There are also other predefined functions allowing you to for example to enhance your presentation with animations of your 3D FES or remove a CV from existing FES. 

Installation:
```bash
pip install metadynminer
```
or
```bash
conda install -c jan8be metadynminer
```

Sample code:

Load your HILLS file: 
```python
hillsfile = metadynminer.Hills(name="HILLS", periodic=[True,True])
```
Compute the free energy surface using the fast Bias Sum Algorithm:
```python
fes = metadynminer.Fes(hillsfile)
```

You can also use slower (but exact) algorithm to sum the hills and compute the free energy surface 
with the option original=True. This algorithm was checked and it gives the same result 
(to the machine level precision) as the PLUMED sum_hills function (for plumed v2.8.0).
```python
fes2 = metadynminer.Fes(hillsfile, original=True)
```

Visualize the free energy surface and save the picture to a file:
```python
fes.plot(png_name="fes.png")
```

Find local minima on the FES, print them and save FES with minima as a picture:
```python
minima = metadynminer.Minima(fes)
print(minima.minima)
minima.plot(png_name="fes.png")
```

You can also plot free energy profile to see, how the differences between each minima were evolving 
during the simulation. Convergence in the free energy profile suggests, that the resulting free energy surface converged to correct values.
```python
fep = metadynminer.FEProfile(minima, hillsfile)
fep.plot(png_name="FEProfile.png")
```

