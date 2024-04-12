# metadynminer.py



[![Build](https://github.com/Jan8be/metadynminer.py/actions/workflows/ci.yml/badge.svg)](https://github.com/Jan8be/metadynminer.py/actions/workflows/ci.yml)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/metadynminer?label=PyPI%20downloads&color=green&link=https%3A%2F%2Fpypi.org%2Fproject%2Fmetadynminer%2F)](https://pypi.org/project/metadynminer/)
[![conda downloads](https://img.shields.io/conda/d/Jan8be/metadynminer?label=Conda%20total%20downloads&color=green&link=https%3A%2F%2Fanaconda.org%2FJan8be%2Fmetadynminer)](https://anaconda.org/Jan8be/metadynminer)


Metadynminer is a package designed to help you analyse output HILLS files from PLUMED metadynamics simulations. 

It is inspired by existing Metadynminer package for R. It supports HILLS files with one, two or three collective variables. 

All built-in functions can be customized with many parameters. You can learn more about that in the [documentation](https://metadynreporter.cz/manual/index.html). 
There are also other predefined functions allowing you to for example to enhance your presentation with animations of your 3D FES or remove a CV from existing FES. 

## Quickstart: run in Binder

Click the icon bellow and wait (couple of minutes) for the container to build and started on public [MyBinder](http://mybinder.org/).

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jan8be/metadynminer.py/main)

Alternatively, for [Metacentrum](https://metacentrum.cz/) users, somewhat better resources are available:

[![Binder](https://binderhub.cloud.e-infra.cz/badge_logo.svg)](https://binderhub.cloud.e-infra.cz/v2/gh/jan8be/metadynminer.py/main?urlpath=lab)

Once in the Jupyterlab environment, upload your ```HILLS``` file and start the ```python_metadynminer.ipynb``` notebook.

## Installation:

```bash
pip install metadynminer
```
or
```bash
conda install -c jan8be metadynminer
```

## Sample code:

Load metadynminer:
```python
import metadynminer
```

Load your HILLS file: 
```python
hillsfile = metadynminer.Hills(name="HILLS")
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

Visualize the free energy surface:
```python
fes.plot()
```

Visualize and save the picture to a file:
```python
fes.plot(png_name="fes.png")
```

Find local minima on the FES, print them and save FES with minima as a picture:
```python
minima = metadynminer.Minima(fes)
print(minima.minima)
minima.plot()
```

You can also plot free energy profile to see, how the differences between each minima were evolving 
during the simulation. Convergence in the free energy profile suggests that the resulting free energy surface converged to correct values.
```python
fep = metadynminer.FEProfile(minima, hillsfile)
fep.plot()
```

