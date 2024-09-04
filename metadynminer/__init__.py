"""
Metadynminer is a package designed to help you analyse output HILLS files from PLUMED metadynamics simulations. 

It is based on Metadynminer package for R programming language, but it is not just a port from R to Python, as it is updated and improved in many aspects. It supports HILLS files with one, two or three collective variables. 

All built-in functions can be easily customized with many parameters. You can learn more about that in the documentation. There are also functions allowing you to enhance your presentation with animations of your 3D FES or remove a CV from existing FES. 
        
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
hillsfile = metadynminer.Hills(name="HILLS")
```

Compute the free energy surface using the fast Bias Sum Algorithm:
```python
fes = metadynminer.Fes(hillsfile)
```

Alternativelly, you can use slower, exact algorithm to sum the hills and compute the free energy surface 
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
fep.plot()
```
"""

name = "metadynminer"
__version__ = "0.8.3"
__author__ = 'Jan Beránek'

__pdoc__ = {}

#print(f"""Welcome to Metadynminer version {__version__}. """)

try:
    import numpy as np
except:
    print("Error while loading numpy")
    exit()
try:
    from matplotlib import pyplot as plt
except:
    print("Error while loading matplotlib pyplot")
    exit()
try:
    from matplotlib import colormaps as cm
except:
    print("Error while loading matplotlib colormaps")
    exit()
try:
    import pandas as pd
except:
    print("Error while loading pandas")
    exit()
try:
    import pyvista as pv
except:
    print("Error while loading pyvista")
    exit()
try:
    import imageio
except:
    print("Error while loading imageio")
    exit()
try:
    import os
except:
    print("Error while loading os")
    exit()
try:
    import copy
except:
    print("Error while loading copy")
    exit()



class TU:
    __pdoc__["TU"] = False
    def __init__(self, name="ps"):
        if name == "ps":
            self.name = name
            self.factor = 1.0
        elif name == "fs":
            self.name = name
            self.factor = 1e-3
        elif name == "ns":
            self.name = name
            self.factor = 1e3
        elif name == "us":
            self.name = "\u03bcs"
            self.factor = 1e6
        elif name == "\u03bcs":
            self.name = "\u03bcs"
            self.factor = 1e6
        elif name == "ms":
            self.name = name
            self.factor = 1e9
        elif name == "s":
            self.name = name
            self.factor = 1e12
        else:
            print(f"Error: invalid name of time unit ({name})")
            
    def inps(self, i):
        if type(i) == int:
            return int(i*self.factor)
        else:
            return i*self.factor

    def intu(self, i):
        if type(i) == int:
            return int(i/self.factor)
        else:
            return i/self.factor



class Hills:
    """
    Object of Hills class are created for loading HILLS files, and obtaining the necessary information from them. 

    Hills files are loaded with command:
    ```python
    hillsfile = metadynminer.Hills(name="HILLS")
    ```
    
    optional parameters:
    
    * name (default="HILLS") = string with name of HILLS file
    
    * ignoretime (default=True) = boolean, if set to False, it will save the time in the HILLS file;
                                if set to True, and timestep is not set, 
                                        each time value will be incremented by the same amount as the time of the first step.
                                        
    * timestep = numeric value of the time difference between hills, in picoseconds
    
    * periodic (default=[False, False]) = list of boolean values telling which CV is periodic.
    
    * cv1per, cv2per, cv3per (defaults = [-numpy.pi, numpy.pi]) = List of two numeric values defining the periodicity of given CV. 
                                        Has to be provided for each periodic CV.

    Hills attributes: 

    * Hills.hills = array of values from the HILLS file
    
    * Hills.cvs = number of CVs

    * Hills.cv1_name = name of CV 1

    * Hills.cv1per = list of two boundaries of the perodicity interval for CV 1

    * Hills.dt = time step between hills in ps

    * Hills.sigma1 = widths of hills along the CV 1 axis

    * Hills.cv1 = list of CV 1 values during simulation

    * Hills.biasf = biasfactors

    * and analogous attributes for other CVs
    """
    
    def __init__(self, name="HILLS", encoding="utf8", ignoretime=True, periodic=None, 
                 cv1per=[-np.pi, np.pi],cv2per=[-np.pi, np.pi],cv3per=[-np.pi, np.pi], timestep=None):
        self.read(name, encoding, ignoretime, periodic, 
                 cv1per=cv1per,cv2per=cv2per,cv3per=cv3per, timestep=timestep)
        self.hillsfilename = name
    
    def read(self, name="HILLS", encoding="utf8", ignoretime=True, periodic=None, 
                 cv1per=[-np.pi, np.pi],cv2per=[-np.pi, np.pi],cv3per=[-np.pi, np.pi], timestep=None):
        with open(name, 'r', encoding=encoding) as hillsfile:
            lines = hillsfile.readlines()
        columns = lines[0].split()
        number_of_columns_head = len(columns) - 2
        
        if number_of_columns_head == 5:
            self.cvs = 1
            self.cv1_name = lines[0].split()[3]
            self.cv1per = cv1per

        elif number_of_columns_head == 7:
            self.cvs = 2
            self.cv1_name = lines[0].split()[3]
            self.cv2_name = lines[0].split()[4]
            self.cv1per = cv1per
            self.cv2per = cv2per
        elif number_of_columns_head == 9:
            self.cvs = 3
            self.cv1_name = lines[0].split()[3]
            self.cv2_name = lines[0].split()[4]
            self.cv3_name = lines[0].split()[5]
            
            self.cv1per = cv1per
            self.cv2per = cv2per
            self.cv3per = cv3per
        else:
            print("Error: Unexpected number of columns in provided HILLS file.")
            return None
        
        self.ignoretime = ignoretime
        if not ignoretime:
            if timestep != None:
                dt = timestep
            else:
                for line in range(len(lines)):
                    if lines[line][0] != "#":
                        dt = round(float(lines[line+1].split()[0]),14) - round(float(lines[line].split()[0]),14)
                        break
        else:
            if timestep != None:
                dt = timestep
            else:
                dt = 1.0
        self.dt = dt
        t = 0
        for line in range(len(lines)):
            if lines[line][0] != "#":
                t += 1 
                if t == 1:
                    if self.cvs == 1:
                        self.sigma1 = float(lines[line].split()[2])
                        self.biasf = float(lines[line].split()[4])
                    elif self.cvs == 2:
                        self.sigma1 = float(lines[line].split()[3])
                        self.sigma2 = float(lines[line].split()[4])
                        self.biasf = float(lines[line].split()[6])
                    elif self.cvs == 3:
                        self.sigma1 = float(lines[line].split()[4])
                        self.sigma2 = float(lines[line].split()[5])
                        self.sigma3 = float(lines[line].split()[6])
                        self.biasf = float(lines[line].split()[8])
                    self.hills = [lines[line].split()]
                    if ignoretime:
                        self.hills[t-1][0] = t*dt
                else:
                    if self.cvs == 1 and len(lines[line].split()) == 5:
                        self.hills.append(lines[line].split())
                        if ignoretime:
                            self.hills[t-1][0] = t*dt
                    if self.cvs == 2 and len(lines[line].split()) == 7:
                        self.hills.append(lines[line].split())
                        if ignoretime:
                            self.hills[t-1][0] = t*dt
                    if self.cvs == 3 and len(lines[line].split()) == 9:
                        self.hills.append(lines[line].split())
                        if ignoretime:
                            self.hills[t-1][0] = t*dt
        
        self.hills = np.array(self.hills, dtype=np.double)
                
        if self.cvs == 1:
            self.cv1 = self.hills[:,1]
            self.sigma1 = self.hills[:,2]
            self.heights = self.hills[:,3]
            self.biasf = self.hills[:,4]
        elif self.cvs == 2:
            self.cv1 = self.hills[:,1]
            self.cv2 = self.hills[:,2]
            self.sigma1 = self.hills[:,3]
            self.sigma2 = self.hills[:,4]
            self.heights = self.hills[:,5]
            self.biasf = self.hills[:,6]
        elif self.cvs == 3:
            self.cv1 = self.hills[:,1]
            self.cv2 = self.hills[:,2]
            self.cv3 = self.hills[:,3]
            self.sigma1 = self.hills[:,4]
            self.sigma2 = self.hills[:,5]
            self.sigma3 = self.hills[:,6]
            self.heights = self.hills[:,7]
            self.biasf = self.hills[:,8]
        print(f"Loaded HILLS file named {name}. ")
        
        # detect periodicity
        if periodic == None:
            if self.cvs >= 1:
                if np.max(np.diff(self.cv1)) > 0.8*(np.max(self.cv1)-np.min(self.cv1)):
                    periodic = list([True])
                else:
                    periodic = list([False])
            if self.cvs >= 2:
                if np.max(np.diff(self.cv2)) > 0.8*(np.max(self.cv2)-np.min(self.cv2)):
                    periodic.append(True)
                else:
                    periodic.append(False)
            if self.cvs >= 3:
                if np.max(np.diff(self.cv3)) > 0.8*(np.max(self.cv3)-np.min(self.cv3)):
                    periodic.append(True)
                else:
                    periodic.append(False)
            print(f"Automatically detected which CVs are periodic: {periodic}. ")
            print("This detection can be overriden by specifying a list of boolean values to 'periodic' keyword. ")
        else:
            if self.cvs == 1:
                if len(periodic) != 1 or (type(periodic[0]) != type(True)):
                    print(f"Error: argument 'periodic' has wrong number of parameters({len(periodic)})")
            if self.cvs == 2:
                if len(periodic) != 2 or (type(periodic[0]) != type(True)) or (type(periodic[1]) != type(True)):
                    print(f"Error: argument 'periodic' has wrong number of parameters({len(periodic)})")
            if self.cvs == 3:
                if len(periodic) != 3 or (type(periodic[0]) != type(True)) or (type(periodic[1]) != type(True)) or (type(periodic[2]) != type(True)):
                    print(f"Error: argument 'periodic' has wrong number of parameters({len(periodic)})")
        
        self.periodic = periodic
        # check periodicity

        if self.cvs >= 1:
            if periodic[0] == False and np.max(np.diff(self.cv1))>0.8*(np.max(self.cv1)-np.min(self.cv1)):
                print(f"WARNING: It looks like CV 1 ({self.cv1_name}) is periodic, however you specified that it is not. This may be completely fine, however, if you forgot to specify the periodicity correctly, it can lead to errors later during FES calculation.")
            if periodic[0] == True and np.max(np.diff(self.cv1))<0.2*(np.max(self.cv1)-np.min(self.cv1)):
                print(f"WARNING: It looks like CV 1 ({self.cv1_name}) is not periodic, however you specified that it is. This may be completely fine, however, if you forgot to specify the periodicity correctly, it can lead to errors later during FES calculation.")
                
        if self.cvs >= 2:
            if periodic[1] == False and np.max(np.diff(self.cv2))>0.8*(np.max(self.cv2)-np.min(self.cv2)):
                print(f"WARNING: It looks like CV 2 ({self.cv2_name}) is periodic, however you specified that it is not. This may be completely fine, however, if you forgot to specify the periodicity correctly, it can lead to errors later during FES calculation.")
            if periodic[1] == True and np.max(np.diff(self.cv2))<0.2*(np.max(self.cv2)-np.min(self.cv2)):
                print(f"WARNING: It looks like CV 2 ({self.cv2_name}) is not periodic, however you specified that it is. This may be completely fine, however, if you forgot to specify the periodicity correctly, it can lead to errors later during FES calculation.")
        
        if self.cvs >= 3:
            if periodic[2] == False and np.max(np.diff(self.cv3))>0.8*(np.max(self.cv3)-np.min(self.cv3)):
                print(f"WARNING: It looks like CV 3 ({self.cv3_name}) is periodic, however you specified that it is not. This may be completely fine, however, if you forgot to specify the periodicity correctly, it can lead to errors later during FES calculation.")
            if periodic[2] == True and np.max(np.diff(self.cv3))<0.2*(np.max(self.cv3)-np.min(self.cv3)):
                print(f"WARNING: It looks like CV 3 ({self.cv3_name}) is not periodic, however you specified that it is. This may be completely fine, however, if you forgot to specify the periodicity correctly, it can lead to errors later during FES calculation.")
                
        return self
    
    
    def get_cv1(self):
        return self.cv1
    def get_cv2(self):
        return self.cv2
    def get_cv3(self):
        return self.cv3
    
    def get_cv1per(self):
        return self.cv1per
    def get_cv2per(self):
        return self.cv2per
    def get_cv3per(self):
        return self.cv3per
    
    def get_periodic(self):
        return self.periodic
    
    def get_cv1_name(self):
        return self.cv1_name
    def get_cv2_name(self):
        return self.cv2_name
    def get_cv3_name(self):
        return self.cv3_name
    
    def get_hills(self):
        return self.hills
    
    def get_number_of_cvs(self):
        return self.cvs
    
    def get_sigma1(self):
        return self.sigma1
    def get_sigma2(self):
        return self.sigma2
    def get_sigma3(self):
        return self.sigma3
    
    def get_heights(self):
        return(self.heights)
    
    __pdoc__["Hills.get_cv1"] = False
    __pdoc__["Hills.get_cv2"] = False
    __pdoc__["Hills.get_cv3"] = False
    __pdoc__["Hills.get_cv1per"] = False
    __pdoc__["Hills.get_cv2per"] = False
    __pdoc__["Hills.get_cv3per"] = False
    __pdoc__["Hills.get_cv1_name"] = False
    __pdoc__["Hills.get_cv2_name"] = False
    __pdoc__["Hills.get_cv3_name"] = False
    __pdoc__["Hills.get_periodic"] = False
    __pdoc__["Hills.get_hills"] = False
    __pdoc__["Hills.get_number_of_cvs"] = False
    __pdoc__["Hills.get_sigma1"] = False
    __pdoc__["Hills.get_sigma2"] = False
    __pdoc__["Hills.get_sigma3"] = False
    __pdoc__["Hills.get_heights"] = False
    __pdoc__["Hills.read"] = False

    def plot_heights(self, png_name=None, energy_unit="kJ/mol", xlabel=None, ylabel=None, label_size=12, image_size=None, image_size_unit="in", dpi=100, tu = "ps", time_min=None, time_max=None, xlim=[None, None], ylim=[None, None], title=None, return_fig=False):
        """
        Function used to visualize heights of the hills that were added during the simulation. 
        
        ```python
        hillsfile.plot_heights(png_name="picture.png")
        ```
        
        Parameters:
        
        * png_name = String. If this parameter is supplied, the picture of FES will be saved under this name to the current working directory.
        
        * energy_unit (default="kJ/mol") = String, used in description of the y axis

        * xlim, ylim (default = [None, None] for both) = list of two values specifying the range of x and y axes respectively. 
        None means that Matplotlib will choose appropriate range, so this keyword is useful when you need to overwrite it. 
        
        * xlabel, ylabel = Strings, if provided, they will be used as labels for the graphs
        
        * labelsize (default = 12) = size of text in labels
        
        * image_size (default = [9,6]) = list of the width and height of the picture

        * image_size_unit (default = "in") = Units for width and height of the picture, accepts "in" as inches, "cm", "mm" or "px" as pixels. 

        * dpi (default = 100) = DPI of the picture

        * tu (default = "ps") = string, time unit to be shown on x axis. Also aplies to time_min and time_max, if those are used. 
        Available options: "s", "ms", "us", "ns", "ps", "fs"

        * time_min, time_max = The time range for plot, closed interval, in the time unit specified by "tu"

        * title = optional, string that defines the title of the graph

        * return_fig (default=False) = whether the method should return the Matplotlib.Pyplot.figure object for further use
        """
        tu = TU(tu)
        
        if time_min == time_max == 0:
            print("Error: Values of start and end time are zero.")
            return None
        if time_min != None:
            #time_min = tu.inps(time_min)
            if tu.inps(time_min) < 0:
                print("Warning: Start time is lower than zero, it will be set to zero instead. ")
                time_min = 0
            if tu.inps(time_min) < int(self.hills[0,0]):
                print(f"Warning: Start time {tu.inps(time_min)} ps is lower than the first time from HILLS file {int(self.hills[0,0])} ps, which will be used instead. ")
                time_min = tu.intu(self.hills[0,0])
        else:
            time_min = tu.intu(self.hills[0,0])
        if time_max != None:
            #time_max = tu.inps(time_max)
            if time_max < time_min:
                print("Warning: End time is lower than start time. Values are flipped. ")
                time_value = time_max
                time_max = time_min
                time_min = time_value
            if tu.inps(time_max) > int(self.hills[-1,0]):
                print(f"Warning: End time {tu.inps(time_max)} ps is higher than number of lines in HILLS file {int(self.hills[-1,0])} ps, which will be used instead. ")
                time_max = tu.intu(self.hills[-1,0])
        else:
            time_max = tu.intu(self.hills[-1,0])
        #print(f"Berofe fes: min {time_min}, max {time_max}")
        
        if not self.ignoretime:
            time_max = int(round(((time_max - self.hills[0,0])/self.dt),0)) + 1
            time_min = int(round(((time_min - self.hills[0,0])/self.dt),0)) + 1

        if image_size == None:
            image_size = [9,6]
        
        if image_size_unit == "cm":
            image_size[0] /= 2.54
            image_size[1] /= 2.54
        elif image_size_unit == "mm":
            image_size[0] /= 25.4
            image_size[1] /= 25.4
        elif image_size_unit == "px":
            image_size[0] /= dpi
            image_size[1] /= dpi
        elif image_size_unit != "in":
            print(f"Warning: unknown image_size_unit value: {image_size_unit}. Using inches instead. ")
            
        fig = plt.figure(figsize=(image_size[0],image_size[1]), dpi=dpi)
        #plt.plot(tu.intu(np.array(range(len(self.heights))))[time_min-1:time_max], self.heights[time_min-1:time_max])
        plt.plot(tu.intu(np.arange(int(self.hills[0,0]), self.heights.shape[0]+1)[int(tu.inps(time_min)-1):int(tu.inps(time_max))]),
                 self.heights[int(tu.inps(time_min)-1):int(tu.inps(time_max))])
        if xlabel == None:
            plt.xlabel(f'time ({tu.name})', size=label_size)
        else:
            plt.xlabel(xlabel, size=label_size)
        if ylabel == None:
            plt.ylabel(f'free energy ({energy_unit})', size=label_size)
        else:
            plt.ylabel(ylabel, size=label_size)
        
        if title != None:
            plt.title(title, size=label_size)

        plt.xticks(fontsize=label_size)
        plt.yticks(fontsize=label_size)
        
        ax = plt.gca()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
            
        if png_name != None:
            plt.savefig(png_name, bbox_inches = 'tight')
        
        if return_fig:
            return fig

    def plot_CV(self, png_name=None, CV=None, xlabel=None, ylabel=None, label_size=12, image_size=None, image_size_unit="in", dpi=100, tu = "ps", time_min=None, time_max=None, points = True, point_size=1, xlim=[None, None], ylim=[None, None], title=None, return_fig=False):
        """
        Function used to visualize CV values from the simulation. 
        
        ```python
        hillsfile.plot_heights(png_name="picture.png")
        ```
        
        Parameters:
        
        * png_name = String. If this parameter is supplied, the picture of FES will be saved under this name to the current working directory.
        
        * CV = Integer, number of the CV to plot (between 1-3)

        * xlim, ylim (default = [None, None] for both) = list of two values specifying the range of x and y axes respectively. 
        None means that Matplotlib will choose appropriate range, so this keyword is useful when you need to overwrite it. 
        
        * xlabel, ylabel = Strings, if provided, they will be used as labels for the graphs
        
        * labelsize (default = 12) = size of text in labels
        
        * image_size (default = [9,6]) = List of the width and height of the picture

        * image_size_unit (default = "in") = Units for width and height of the picture, accepts "in" as inches, "cm", "mm" or "px" as pixels. 

        * dpi (default = 100) = DPI of the resulting image. 
        
        * tu (default = "ps") = string, time unit to be shown on x axis. Also aplies to time_min and time_max, if those are used. Available options: "s", "ms", "us", "ns", "ps", "fs"

        * time_min, time_max = The time range for plot

        * points (default=True) = Boolean value; if True, plot type will be scatter plot, which is better for periodic CVs; if False, it will be line plot, which is sometimes more suitable for non-periodic CVs. 

        * point_size (default = 1) = The size of dots in the plot

        * title = optional, string that defines the title of the graph

        * return_fig (default=False) = whether the method should return the Matplotlib.Pyplot.figure object for further use
        """
        if CV==None:
            print("Error: CV was not chosen")
            return None
        if CV>self.cvs:
            print(f"Error: CV {CV} is not available. ")
            return None
        if CV==1.0:
            CV=1
        if CV==2.0:
            CV=2
        if CV==3.0:
            CV=3
        if CV!=1 and CV!=2 and CV!=3:
            print(f"Error: supplied value of CV {CV} is not correct value")
            return None
            
        tu = TU(tu)
        
        if time_min == time_max == 0:
            print("Error: Values of start and end time are zero.")
            return None
        if time_min != None:
            #time_min = tu.inps(time_min)
            if tu.inps(time_min) < 0:
                print("Warning: Start time is lower than zero, it will be set to zero instead. ")
                time_min = 0
            if tu.inps(time_min) < int(self.hills[0,0]):
                print(f"Warning: Start time {tu.inps(time_min)} ps is lower than the first time from HILLS file {int(self.hills[0,0])} ps, which will be used instead. ")
                time_min = tu.intu(self.hills[0,0])
        else:
            time_min = tu.intu(self.hills[0,0])
        if time_max != None:
            #time_max = tu.inps(time_max)
            if time_max < time_min:
                print("Warning: End time is lower than start time. Values are flipped. ")
                time_value = time_max
                time_max = time_min
                time_min = time_value
            if tu.inps(time_max) > int(self.hills[-1,0]):
                print(f"Warning: End time {tu.inps(time_max)} ps is higher than number of lines in HILLS file {int(self.hills[-1,0])} ps, which will be used instead. ")
                time_max = tu.intu(self.hills[-1,0])
        else:
            time_max = tu.intu(self.hills[-1,0])
        
        if not self.ignoretime:
            time_max = int(round(((time_max - self.hills[0,0])/self.dt),0)) + 1
            time_min = int(round(((time_min - self.hills[0,0])/self.dt),0)) + 1

        if image_size == None:
            image_size = [9,6]
        
        if image_size_unit == "cm":
            image_size[0] /= 2.54
            image_size[1] /= 2.54
        elif image_size_unit == "mm":
            image_size[0] /= 25.4
            image_size[1] /= 25.4
        elif image_size_unit == "px":
            image_size[0] /= dpi
            image_size[1] /= dpi
        elif image_size_unit != "in":
            print(f"Warning: unknown image_size_unit value: {image_size_unit}. Using inches instead. ")
        
        fig = plt.figure(figsize=(image_size[0],image_size[1]), dpi=dpi)
        if points:
            if CV==1:
                plt.scatter(tu.intu(np.array(range(int(round(tu.inps(time_min),0)),int(round(tu.inps(time_max)+1,0)),int(round(self.dt,0))))), 
                            self.cv1[int(tu.inps(time_min)-1):int(tu.inps(time_max))], s=point_size)
            if CV==2:
                plt.scatter(tu.intu(np.array(range(int(round(tu.inps(time_min),0)),int(round(tu.inps(time_max)+1,0)),int(round(self.dt,0))))), 
                            self.cv2[int(tu.inps(time_min)-1):int(tu.inps(time_max))], s=point_size)
            if CV==3:
                plt.scatter(tu.intu(np.array(range(int(round(tu.inps(time_min),0)),int(round(tu.inps(time_max)+1,0)),int(round(self.dt,0))))), 
                            self.cv3[int(tu.inps(time_min)-1):int(tu.inps(time_max))], s=point_size)
        else:
            if CV==1:
                plt.plot(tu.intu(np.array(range(int(round(tu.inps(time_min),0)),int(round(tu.inps(time_max)+1,0)),int(round(self.dt,0))))), 
                         self.cv1[int(tu.inps(time_min)-1):int(tu.inps(time_max))], lw=point_size)
            if CV==2:
                plt.plot(tu.intu(np.array(range(int(round(tu.inps(time_min),0)),int(round(tu.inps(time_max)+1,0)),int(round(self.dt,0))))), 
                         self.cv2[int(tu.inps(time_min)-1):int(tu.inps(time_max))], lw=point_size)
            if CV==3:
                plt.plot(tu.intu(np.array(range(int(round(tu.inps(time_min),0)),int(round(tu.inps(time_max)+1,0)),int(round(self.dt,0))))), 
                         self.cv3[int(tu.inps(time_min)-1):int(tu.inps(time_max))], lw=point_size)
            
        if xlabel == None:
            plt.xlabel(f'time ({tu.name})', size=label_size)
        else:
            plt.xlabel(xlabel, size=label_size)
        if ylabel == None:
            if CV==1:
                plt.ylabel(f'CV {CV} - {self.cv1_name}', size=label_size)
            if CV==2:
                plt.ylabel(f'CV {CV} - {self.cv2_name}', size=label_size)
            if CV==3:
                plt.ylabel(f'CV {CV} - {self.cv3_name}', size=label_size)
        else:
            plt.ylabel(ylabel, size=label_size)
        
        if title != None:
            plt.title(title, size=label_size)

        plt.xticks(fontsize=label_size)
        plt.yticks(fontsize=label_size)
                
        ax = plt.gca()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
            
        if png_name != None:
            plt.savefig(png_name, bbox_inches = 'tight')
        
        if return_fig:
            return fig

class Fes: 
    """
    Object of this class is created to compute the free energy surface corresponding to the provided Hills object. 
    Command:
    ```python
    fes = metadynminer.Fes(hills=hillsfile)
    ```
    parameters:
    
    * hills = Hills object
    
    * resolution (default=256) = should be positive integer, controls the resolution of FES
    
    * original (default=False) = boolean, if False, FES will be calculated using very fast, but not
    'exact' Bias Sum Algorithm
                                        if True, FES will be calculated with slower algorithm, but it will be exactly the same as FES calculated 
                                        with PLUMED sum_hills function
                                        
    * cv1range, cv2range, cv3range = lists of two numbers, defining lower and upper bound of the respective CV (in the units of the CVs)
    
    * time_min, time_max = Limit points for closed interval of times from the HILLS file from which the FES will be constructed. Useful for making animations of flooding of the FES during simulation. The values here should have the same units as those in the HILLS file, especially if ignoretime = False. 

    * tu (default = "ps") = string, time unit for time_min and time_max parameters, if those are used. Available options: "s", "ms", "us", "ns", "ps", "fs"

    * print_output (default = True), tells whether the function should print messages like progress of calculation

    * subtract_min (deault = True), whether the global minimum value should be subtracted from the whole FES, so that global minimum has zero free energy and other states have higher values

    Fes atributes: 

    * Fes.fes = free energy surface
    
    * Fes.res = resolution of the FES

    * Fes.original = whether FES should be calculated by plumed's original precise algorithm

    * Fes.cvs = number of CVs

    * Fes.hills = array with values from the HILLS file

    * Fes.periodic = list of boolean values setting which CVs are periodic

    * Fes.cv1min, Fes.cv1max = minimum and maximum values on the FES, respectively

    * Fes.cv1_fes_range, Fes.cv2_fes_range, Fes.cv3_fes_range = FES range = maximum - minimum values on FES for given CV

    * Fes.cv1_name, Fes.cv2_name, Fes.cv3_name = CV names
    
    """
    
    __pdoc__["Fes.makefes"] = False
    __pdoc__["Fes.makefes2"] = False
    __pdoc__["Fes.set_fes"] = False
    
    def __init__(self, hills=None, resolution=256, original=False, \
                 calculate_new_fes=True, cv1range=None, cv2range=None, cv3range=None, \
                 time_min=None, time_max=None, tu="ps", print_output=True, subtract_min = True):
        self.res = resolution
        self.original = original
        self.cv1range = cv1range
        self.cv2range = cv2range
        self.cv3range = cv3range
        if hills != None:
            self.hills = hills
            self.cvs = hills.get_number_of_cvs()
            self.heights = hills.get_heights()
            self.periodic = hills.get_periodic()
            self.biasf = hills.biasf
            self.ignoretime = hills.ignoretime
            self.dt = hills.dt
            
            if cv1range!=None and len(cv1range) != 2:
                print("Error: You have to specify CV1 range as a list of two values. ")
            if cv2range!=None and len(cv2range) != 2:
                print("Error: You have to specify CV2 range as a list of two values. ")
            if cv3range!=None and len(cv3range) != 2:
                print("Error: You have to specify CV3 range as a list of two values. ")
            
            if self.cvs >= 1:
                self.cv1 = hills.get_cv1()
                self.s1 = hills.get_sigma1()

                self.cv1_name = hills.get_cv1_name()
                self.cv1per = hills.get_cv1per()

                self.fes = np.zeros((resolution))

                if self.periodic[0]:
                    cv1min = self.cv1per[0]
                    cv1max = self.cv1per[1]
                    if cv1range != None:
                        print("Warning: CV1 is specified as periodic, to change it's range you must specify the periodicity to cv1per parameter when loading the HILLS file. The cv1range parameter is ignored. ")
                else:
                    if cv1range == None:
                        cv1min = np.min(self.cv1) - 1e-8
                        cv1max = np.max(self.cv1) + 1e-8
                        cv1min -= (cv1max-cv1min)*0.15          
                        cv1max += (cv1max-cv1min)*0.15
                    else:
                        if cv1range[1] <= cv1range[0]:
                            print(f"Error: Wrong values of cv1range: {cv1range}. ")
                            return None
                        cv1min = cv1range[0]
                        cv1max = cv1range[1]
                self.cv1min = cv1min
                self.cv1max = cv1max
                self.cv1_fes_range = self.cv1max - self.cv1min
                
                
                if not original:
                    if ((np.max(self.s1)/np.min(self.s1))>1.00000001):
                        print("""Error: Bias sum algorithm only works for hills files 
                        in which all hills have the same width. 
                        For this file, you need the slower but exact, algorithm, to use it, 
                        specify the argument 'original=True'.""")
                        return None

            if self.cvs >= 2:
                self.cv2 = hills.get_cv2()
                self.s2 = hills.get_sigma2()
                
                self.cv2_name = hills.get_cv2_name()
                self.cv2per = hills.get_cv2per()

                self.fes = np.zeros((resolution, resolution))

                if self.periodic[1]:
                    cv2min = self.cv2per[0]
                    cv2max = self.cv2per[1]
                    if cv2range != None:
                        print("Warning: CV2 is specified as periodic, to change it's range you must specify the periodicity to cv2per parameter when loading the HILLS file. The cv2range parameter is ignored. ")
                else:
                    if cv2range == None:
                        cv2min = np.min(self.cv2) - 1e-8
                        cv2max = np.max(self.cv2) + 1e-8
                        cv2min -= (cv2max-cv2min)*0.15          
                        cv2max += (cv2max-cv2min)*0.15
                    else:
                        if cv2range[1] <= cv2range[0]:
                            print(f"Error: Wrong values of cv2range: {cv2range}. ")
                            return None
                        cv2min = cv2range[0]
                        cv2max = cv2range[1]
                self.cv2min = cv2min
                self.cv2max = cv2max
                self.cv2_fes_range = self.cv2max - self.cv2min


                
                if not original:
                    if ((np.max(self.s2)/np.min(self.s2))>1.00000001):
                        print("""Error: Bias sum algorithm only works for hills files 
                        in which all hills have the same width. 
                        For this file, you need the slower, but exact, algorithm, 
                        specify the argument 'original=True'.""")
                        return None

            if self.cvs == 3:
                self.cv3 = hills.get_cv3()
                self.s3 = hills.get_sigma3()

                self.cv3_name = hills.get_cv3_name()
                self.cv3per = hills.get_cv3per()

                self.fes = np.zeros((resolution, resolution, resolution))
                
                if self.periodic[2]:
                    cv3min = self.cv3per[0]
                    cv3max = self.cv3per[1]
                    if cv3range != None:
                        print("Warning: CV3 is specified as periodic, to change it's range you must specify the periodicity to cv3per parameter when loading the HILLS file. The cv3range is ignored. ")
                else:
                    if cv3range == None:
                        cv3min = np.min(self.cv3) - 1e-8
                        cv3max = np.max(self.cv3) + 1e-8
                        cv3min -= (cv3max-cv3min)*0.15          
                        cv3max += (cv3max-cv3min)*0.15
                    else:
                        if cv3range[1] <= cv3range[0]:
                            print(f"Error: Wrong values of cv3range: {cv3range}. ")
                            return None
                        cv3min = cv3range[0]
                        cv3max = cv3range[1]
                self.cv3min = cv3min
                self.cv3max = cv3max
                self.cv3_fes_range = self.cv3max - self.cv3min
                
                if not original:
                    if ((np.max(self.s3)/np.min(self.s3))>1.00000001):
                        print("""Error: Bias sum algorithm only works for hills files 
                        in which all hills have the same width of given CV. 
                        For this file, you need the exact algorithm, to do that, 
                        specify the argument 'original=True'.""")
                        return None
                        
            tu = TU(tu)
        
            if time_min == time_max == 0:
                print("Error: Values of start and end time are zero.")
                return None
            if time_min != None:
                #time_min = tu.inps(time_min)
                if time_min < 0:
                    print("Warning: Start time is lower than zero, it will be set to zero instead. ")
                    time_min = 0
                if tu.inps(time_min) < int(hills.hills[0,0]):
                    print(f"Warning: Start time {tu.inps(time_min)} is lower than the first time from HILLS file {int(hills.hills[0,0])}, which will be used instead. ")
                    time_min = tu.intu(int(hills.hills[0,0]))
            else:
                time_min = tu.intu(int(hills.hills[0,0]))
            if time_max != None:
                #time_max = int(tu.inps(time_max))
                if time_max < time_min:
                    print("Warning: End time is lower than start time. Values are flipped. ")
                    time_value = time_max
                    time_max = time_min
                    time_min = time_value
                if tu.inps(time_max) > int(hills.hills[-1,0]):
                    print(f"Warning: End time {tu.inps(time_max)} is higher than number of lines in HILLS file {int(hills.hills[-1,0])}, which will be used instead. ")
                    time_max = tu.intu(int(hills.hills[-1,0]))
            else:
                time_max = tu.intu(int(hills.hills[-1,0]))
            
            if not self.ignoretime:
                time_max = int(round(((time_max - hills.hills[0,0])/self.dt),0)) + 1
                time_min = int(round(((time_min - hills.hills[0,0])/self.dt),0)) + 1
                #print(f"Berofe fes: min {time_min}, max {time_max}")

            if calculate_new_fes:
                if original:
                    self.makefes2(resolution, int(tu.inps(time_min)), int(tu.inps(time_max)), 
                                  print_output=print_output, subtract_min = subtract_min)
                else:
                    self.makefes(resolution, int(tu.inps(time_min)), int(tu.inps(time_max)), 
                                 print_output=print_output, subtract_min = subtract_min)
        
    def makefes(self, resolution, time_min, time_max, print_output=True, subtract_min = True):
        """
        Function used internally for summing hills in Hills object with the fast Bias Sum Algorithm. 
        """
        #self.res = resolution
        #if self.res % 2 == 0:
        #    self.res += 1
        #print(f"min: {time_min}, max: {time_max}")
        
        if self.cvs == 1:
            cv1bin = np.ceil((self.cv1-self.cv1min)*self.res/(self.cv1_fes_range))
            cv1bin = cv1bin.astype(int)
            
            s1res = (self.s1[0]*self.res)/(self.cv1_fes_range)
            
            self.cv1bin = cv1bin
            
            gauss_res = 8*s1res
            gauss_res = int(gauss_res)
            if gauss_res%2 == 0:
                gauss_res += 1
            
            gauss_center = int((gauss_res-1)/2)+1
            gauss = np.zeros((gauss_res))
            for i in range(gauss_res):
                gauss[int(i)] = -np.exp(-((i+1)-gauss_center)**2/(2*s1res**2))
                
            fes = np.zeros((self.res))
            
            for line in range(time_min-1, time_max):
                if print_output and (((line) % 500 == 0) or (line == time_max-1)):
                    print(f"Constructing free energy surface: {((line+1-time_min+1)/(time_max-time_min+1)):.1%} finished", end="\r")
                
                gauss_center_to_end = int((gauss_res-1)/2)
                
                fes_to_edit_cv1 = [cv1bin[line]-1-gauss_center_to_end,
                                   cv1bin[line]-1+gauss_center_to_end]
                fes_crop_cv1 = [max(0,fes_to_edit_cv1[0]),min(self.res-1,fes_to_edit_cv1[1])]
                
                gauss_crop_cv1 = [max(0,gauss_center_to_end-(cv1bin[line]-1)),
                                      gauss_res-1-max(0,(cv1bin[line]-1)+gauss_center_to_end-self.res+1)]
                
                fes[fes_crop_cv1[0]:fes_crop_cv1[1]+1]\
                        += gauss[gauss_crop_cv1[0]:gauss_crop_cv1[1]+1]\
                        * self.heights[line]
                
                if self.periodic[0]:
                    if cv1bin[line] < gauss_center_to_end:
                        fes_crop_cv1_p = [self.res-1+(cv1bin[line]-gauss_center_to_end),self.res-1]
                        gauss_crop_cv1_p = [0,gauss_center_to_end-cv1bin[line]]
                        fes[fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1]\
                                    += gauss[gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1]\
                                    * self.heights[line]
                    
                    if cv1bin[line] > (self.res-gauss_center_to_end):
                        fes_crop_cv1_p = [0,gauss_center_to_end+cv1bin[line]-self.res-1]
                        gauss_crop_cv1_p = [gauss_res-(gauss_center_to_end+cv1bin[line]-self.res),gauss_res-1]
                        fes[fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1]\
                                    += gauss[gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1]\
                                    * self.heights[line]
            if print_output:
                print("\n")
            if subtract_min:
                fes = fes-np.min(fes)
            self.fes = np.array(fes)
            
        elif self.cvs == 2:
            cv1bin = np.ceil((self.cv1-self.cv1min)*self.res/(self.cv1_fes_range))
            cv2bin = np.ceil((self.cv2-self.cv2min)*self.res/(self.cv2_fes_range))
                        
            cv1bin = cv1bin.astype(int)
            cv2bin = cv2bin.astype(int)
            
            s1res = (self.s1[0]*self.res)/(self.cv1_fes_range)
            s2res = (self.s2[0]*self.res)/(self.cv2_fes_range)
            
            gauss_res = max(8*s1res, 8*s2res)
            gauss_res = int(gauss_res)
            if gauss_res%2 == 0:
                gauss_res += 1
            
            gauss_center = int((gauss_res-1)/2)+1
            gauss = np.zeros((gauss_res,gauss_res))
            for i in range(gauss_res):
                for j in range(gauss_res):
                    #dp2 = ((i+1)-gauss_center)**2/(2*s1res**2) + ((j+1)-gauss_center)**2/(2*s2res**2)
                    #if dp2 < 6.25:
                    #    gauss[int(i), int(j)] = -np.exp(-dp2) * 1.00193418799744762399 - 0.00193418799744762399
                    gauss[int(i), int(j)] = -np.exp(-(((i+1)-gauss_center)**2/(2*s1res**2) + ((j+1)-gauss_center)**2/(2*s2res**2)))
            
            fes = np.zeros((self.res,self.res))
            for line in range(time_min-1, time_max):
                if print_output and (((line) % 500 == 0) or (line == time_max-1)):
                    print(f"Constructing free energy surface: {((line+1-time_min+1)/(time_max-time_min+1)):.1%} finished.", end="\r")
                
                #fes_center = int((self.res-1)/2)
                gauss_center_to_end = int((gauss_res-1)/2)
                #print(f"\ng_res: {gauss_res}, gauss_center_to_end: {gauss_center_to_end}, cvib:{cv1bin[line]}, cv2b:{cv2bin[line]}")
                
                fes_to_edit_cv1 = [cv1bin[line]-1-gauss_center_to_end,
                                   cv1bin[line]-1+gauss_center_to_end]
                fes_to_edit_cv2 = [cv2bin[line]-1-gauss_center_to_end,
                                   cv2bin[line]-1+gauss_center_to_end]
                fes_crop_cv1 = [max(0,fes_to_edit_cv1[0]),min(self.res-1,fes_to_edit_cv1[1])]
                fes_crop_cv2 = [max(0,fes_to_edit_cv2[0]),min(self.res-1,fes_to_edit_cv2[1])]
                
                gauss_crop_cv1 = [max(0,gauss_center_to_end-(cv1bin[line]-1)),
                                      gauss_res-1-max(0,(cv1bin[line]-1)+gauss_center_to_end-self.res+1)]
                gauss_crop_cv2 = [max(0,gauss_center_to_end-(cv2bin[line]-1)),
                                      gauss_res-1-max(0,(cv2bin[line]-1)+gauss_center_to_end-self.res+1)]
                
                fes[fes_crop_cv1[0]:fes_crop_cv1[1]+1,fes_crop_cv2[0]:fes_crop_cv2[1]+1]\
                        += gauss[gauss_crop_cv1[0]:gauss_crop_cv1[1]+1,gauss_crop_cv2[0]:gauss_crop_cv2[1]+1]\
                        * self.heights[line]
                
                if self.periodic[0]:
                    if cv1bin[line] < gauss_center_to_end:
                        fes_crop_cv1_p = [self.res-1+(cv1bin[line]-gauss_center_to_end),self.res-1]
                        gauss_crop_cv1_p = [0,gauss_center_to_end-cv1bin[line]]
                        fes[fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1,fes_crop_cv2[0]:fes_crop_cv2[1]+1]\
                                    += gauss[gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1,gauss_crop_cv2[0]:gauss_crop_cv2[1]+1]\
                                    * self.heights[line]
                    
                    if cv1bin[line] > (self.res-gauss_center_to_end):
                        fes_crop_cv1_p = [0,gauss_center_to_end+cv1bin[line]-self.res-1]
                        gauss_crop_cv1_p = [gauss_res-(gauss_center_to_end+cv1bin[line]-self.res),gauss_res-1]
                        fes[fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1,fes_crop_cv2[0]:fes_crop_cv2[1]+1]\
                                    += gauss[gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1,gauss_crop_cv2[0]:gauss_crop_cv2[1]+1]\
                                    * self.heights[line]
                
                if self.periodic[1]:
                    if cv2bin[line] < gauss_center_to_end:
                        fes_crop_cv2_p = [self.res-1+(cv2bin[line]-gauss_center_to_end),self.res-1]
                        gauss_crop_cv2_p = [0,gauss_center_to_end-cv2bin[line]]
                        fes[fes_crop_cv1[0]:fes_crop_cv1[1]+1,fes_crop_cv2_p[0]:fes_crop_cv2_p[1]+1]\
                                    += gauss[gauss_crop_cv1[0]:gauss_crop_cv1[1]+1,gauss_crop_cv2_p[0]:gauss_crop_cv2_p[1]+1]\
                                    * self.heights[line]
                    
                    if cv2bin[line] > (self.res-gauss_center_to_end):
                        fes_crop_cv2_p = [0,gauss_center_to_end+cv2bin[line]-self.res-1]
                        gauss_crop_cv2_p = [gauss_res-(gauss_center_to_end+cv2bin[line]-self.res),gauss_res-1]
                        fes[fes_crop_cv1[0]:fes_crop_cv1[1]+1,fes_crop_cv2_p[0]:fes_crop_cv2_p[1]+1]\
                                    += gauss[gauss_crop_cv1[0]:gauss_crop_cv1[1]+1,gauss_crop_cv2_p[0]:gauss_crop_cv2_p[1]+1]\
                                    * self.heights[line]
                
                if self.periodic[0] and self.periodic[1]:
                    if ((cv1bin[line] < gauss_center_to_end) or (cv1bin[line] > (self.res-gauss_center_to_end))) \
                            and ((cv2bin[line] < gauss_center_to_end) or (cv2bin[line] > (self.res-gauss_center_to_end))):
                        fes[fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1,fes_crop_cv2_p[0]:fes_crop_cv2_p[1]+1]\
                                    += gauss[gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1,gauss_crop_cv2_p[0]:gauss_crop_cv2_p[1]+1]\
                                    * self.heights[line]
            if print_output:
                print("\n")
            if subtract_min:
                fes = fes-np.min(fes)
            self.fes = np.array(fes)
        elif self.cvs == 3:
            cv1bin = np.ceil((self.cv1-self.cv1min)*self.res/(self.cv1_fes_range))
            cv2bin = np.ceil((self.cv2-self.cv2min)*self.res/(self.cv2_fes_range))
            cv3bin = np.ceil((self.cv3-self.cv3min)*self.res/(self.cv3_fes_range))
                        
            cv1bin = cv1bin.astype(int)
            cv2bin = cv2bin.astype(int)
            cv3bin = cv3bin.astype(int)
            
            s1res = (self.s1[0]*self.res)/(self.cv1_fes_range)
            s2res = (self.s2[0]*self.res)/(self.cv2_fes_range)
            s3res = (self.s3[0]*self.res)/(self.cv3_fes_range)
            
            gauss_res = max(10*s1res, 10*s2res, 10*s3res)
            gauss_res = int(gauss_res)
            if gauss_res%2 == 0:
                gauss_res += 1
            gauss_center = int((gauss_res-1)/2)+1
            gauss = np.zeros((gauss_res,gauss_res,gauss_res))
            for i in range(gauss_res):
                for j in range(gauss_res):
                    for k in range(gauss_res):
                        #dp2 = ((i+1)-gauss_center)**2/(2*s1res**2) + ((j+1)-gauss_center)**2/(2*s2res**2) + ((k+1)-gauss_center)**2/(2*s3res**2)
                        #if dp2 < 6.25:
                        #    gauss[int(i), int(j), int(k)] = -np.exp(-dp2) * 1.00193418799744762399 - 0.00193418799744762399
                        gauss[int(i), int(j), int(k)] = -np.exp(-(((i+1)-gauss_center)**2/(2*s1res**2) + 
                                                                  ((j+1)-gauss_center)**2/(2*s2res**2) + 
                                                                  ((k+1)-gauss_center)**2/(2*s3res**2)))
            
            fes = np.zeros((self.res, self.res, self.res))
            
            for line in range(time_min-1, time_max):
                if print_output and (((line) % 500 == 0) or (line == time_max-1)):
                    print(f"Constructing free energy surface: {((line+1-time_min+1)/(time_max-time_min+1)):.1%} finished", end="\r")
                
                #fes_center = int((self.res-1)/2)
                gauss_center_to_end = int((gauss_res-1)/2)
                
                fes_to_edit_cv1 = [cv1bin[line]-1-gauss_center_to_end,
                                   cv1bin[line]-1+gauss_center_to_end]
                fes_to_edit_cv2 = [(cv2bin[line]-1)-gauss_center_to_end,
                                   (cv2bin[line]-1)+gauss_center_to_end]
                fes_to_edit_cv3 = [(cv3bin[line]-1)-gauss_center_to_end,
                                   (cv3bin[line]-1)+gauss_center_to_end]
                fes_crop_cv1 = [max(0,fes_to_edit_cv1[0]),min(self.res-1,fes_to_edit_cv1[1])]
                fes_crop_cv2 = [max(0,fes_to_edit_cv2[0]),min(self.res-1,fes_to_edit_cv2[1])]
                fes_crop_cv3 = [max(0,fes_to_edit_cv3[0]),min(self.res-1,fes_to_edit_cv3[1])]
                
                gauss_crop_cv1 = [max(0,gauss_center_to_end-(cv1bin[line]-1)),
                                      gauss_res-1-max(0,(cv1bin[line]-1)+gauss_center_to_end-self.res+1)]
                gauss_crop_cv2 = [max(0,gauss_center_to_end-(cv2bin[line]-1)),
                                      gauss_res-1-max(0,(cv2bin[line]-1)+gauss_center_to_end-self.res+1)]
                gauss_crop_cv3 = [max(0,gauss_center_to_end-(cv3bin[line]-1)),
                                      gauss_res-1-max(0,(cv3bin[line]-1)+gauss_center_to_end-self.res+1)]
                
                fes[fes_crop_cv1[0]:fes_crop_cv1[1]+1,\
                    fes_crop_cv2[0]:fes_crop_cv2[1]+1,\
                    fes_crop_cv3[0]:fes_crop_cv3[1]+1]\
                        += gauss[gauss_crop_cv1[0]:gauss_crop_cv1[1]+1,\
                                 gauss_crop_cv2[0]:gauss_crop_cv2[1]+1,\
                                 gauss_crop_cv3[0]:gauss_crop_cv3[1]+1]\
                        * self.heights[line]
                
                if self.periodic[0]:
                    if cv1bin[line] < gauss_center_to_end:
                        fes_crop_cv1_p = [self.res-1+(cv1bin[line]-gauss_center_to_end),self.res-1]
                        gauss_crop_cv1_p = [0,gauss_center_to_end-cv1bin[line]]
                        fes[fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1,\
                            fes_crop_cv2[0]:fes_crop_cv2[1]+1,\
                            fes_crop_cv3[0]:fes_crop_cv3[1]+1]\
                                    += gauss[gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1,\
                                             gauss_crop_cv2[0]:gauss_crop_cv2[1]+1,\
                                             gauss_crop_cv3[0]:gauss_crop_cv3[1]+1]\
                                    * self.heights[line]
                    
                    if cv1bin[line] > (self.res-gauss_center_to_end):
                        fes_crop_cv1_p = [0,gauss_center_to_end+cv1bin[line]-self.res-1]
                        gauss_crop_cv1_p = [gauss_res-(gauss_center_to_end+cv1bin[line]-self.res),gauss_res-1]
                        fes[fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1,\
                            fes_crop_cv2[0]:fes_crop_cv2[1]+1,\
                            fes_crop_cv3[0]:fes_crop_cv3[1]+1]\
                                    += gauss[gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1,\
                                             gauss_crop_cv2[0]:gauss_crop_cv2[1]+1,\
                                             gauss_crop_cv3[0]:gauss_crop_cv3[1]+1]\
                                    * self.heights[line]
                
                if self.periodic[1]:
                    if cv2bin[line] < gauss_center_to_end:
                        fes_crop_cv2_p = [self.res-1+(cv2bin[line]-gauss_center_to_end),self.res-1]
                        gauss_crop_cv2_p = [0,gauss_center_to_end-cv2bin[line]]
                        fes[fes_crop_cv1[0]:fes_crop_cv1[1]+1,\
                            fes_crop_cv2_p[0]:fes_crop_cv2_p[1]+1,\
                            fes_crop_cv3[0]:fes_crop_cv3[1]+1]\
                                    += gauss[gauss_crop_cv1[0]:gauss_crop_cv1[1]+1,\
                                             gauss_crop_cv2_p[0]:gauss_crop_cv2_p[1]+1,\
                                             gauss_crop_cv3[0]:gauss_crop_cv3[1]+1]\
                                    * self.heights[line]
                    
                    if cv2bin[line] > (self.res-gauss_center_to_end):
                        fes_crop_cv2_p = [0,gauss_center_to_end+cv2bin[line]-self.res-1]
                        gauss_crop_cv2_p = [gauss_res-(gauss_center_to_end+cv2bin[line]-self.res),gauss_res-1]
                        fes[fes_crop_cv1[0]:fes_crop_cv1[1]+1,\
                            fes_crop_cv2_p[0]:fes_crop_cv2_p[1]+1,\
                            fes_crop_cv3[0]:fes_crop_cv3[1]+1]\
                                    += gauss[gauss_crop_cv1[0]:gauss_crop_cv1[1]+1,\
                                             gauss_crop_cv2_p[0]:gauss_crop_cv2_p[1]+1,\
                                             gauss_crop_cv3[0]:gauss_crop_cv3[1]+1]\
                                    * self.heights[line]
                if self.periodic[2]:
                    if cv3bin[line] < gauss_center_to_end:
                        fes_crop_cv3_p = [self.res-1+(cv3bin[line]-gauss_center_to_end),self.res-1]
                        gauss_crop_cv3_p = [0,gauss_center_to_end-cv3bin[line]]
                        fes[fes_crop_cv1[0]:fes_crop_cv1[1]+1,\
                            fes_crop_cv2[0]:fes_crop_cv2[1]+1,\
                            fes_crop_cv3_p[0]:fes_crop_cv3_p[1]+1]\
                                    += gauss[gauss_crop_cv1[0]:gauss_crop_cv1[1]+1,\
                                             gauss_crop_cv2[0]:gauss_crop_cv2[1]+1,\
                                             gauss_crop_cv3_p[0]:gauss_crop_cv3_p[1]+1]\
                                    * self.heights[line]
                    
                    if cv3bin[line] > (self.res-gauss_center_to_end):
                        fes_crop_cv3_p = [0,gauss_center_to_end+cv3bin[line]-self.res-1]
                        gauss_crop_cv3_p = [gauss_res-(gauss_center_to_end+cv3bin[line]-self.res),gauss_res-1]
                        fes[fes_crop_cv1[0]:fes_crop_cv1[1]+1,\
                            fes_crop_cv2[0]:fes_crop_cv2[1]+1,\
                            fes_crop_cv3_p[0]:fes_crop_cv3_p[1]+1]\
                                    += gauss[gauss_crop_cv1[0]:gauss_crop_cv1[1]+1,\
                                             gauss_crop_cv2[0]:gauss_crop_cv2[1]+1,\
                                             gauss_crop_cv3_p[0]:gauss_crop_cv3_p[1]+1]\
                                    * self.heights[line]
                
                if self.periodic[0] and self.periodic[1]:
                    if ((cv1bin[line] < gauss_center_to_end) or (cv1bin[line] > (self.res-gauss_center_to_end))) \
                            and ((cv2bin[line] < gauss_center_to_end) or (cv2bin[line] > (self.res-gauss_center_to_end))):
                        fes[fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1,\
                            fes_crop_cv2_p[0]:fes_crop_cv2_p[1]+1,\
                            fes_crop_cv3[0]:fes_crop_cv3[1]+1]\
                                    += gauss[gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1,\
                                             gauss_crop_cv2_p[0]:gauss_crop_cv2_p[1]+1,\
                                             gauss_crop_cv3[0]:gauss_crop_cv3[1]+1]\
                                    * self.heights[line]
                
                if self.periodic[0] and self.periodic[2]:
                    if ((cv1bin[line] < gauss_center_to_end) or (cv1bin[line] > (self.res-gauss_center_to_end))) \
                            and ((cv3bin[line] < gauss_center_to_end) or (cv3bin[line] > (self.res-gauss_center_to_end))):
                        fes[fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1,\
                            fes_crop_cv2[0]:fes_crop_cv2[1]+1,\
                            fes_crop_cv3_p[0]:fes_crop_cv3_p[1]+1]\
                                    += gauss[gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1,\
                                             gauss_crop_cv2[0]:gauss_crop_cv2[1]+1,\
                                             gauss_crop_cv3_p[0]:gauss_crop_cv3_p[1]+1]\
                                    * self.heights[line]
                
                if self.periodic[1] and self.periodic[2]:
                    if ((cv2bin[line] < gauss_center_to_end) or (cv2bin[line] > (self.res-gauss_center_to_end))) \
                            and ((cv3bin[line] < gauss_center_to_end) or (cv3bin[line] > (self.res-gauss_center_to_end))):
                        fes[fes_crop_cv1[0]:fes_crop_cv1[1]+1,\
                            fes_crop_cv2_p[0]:fes_crop_cv2_p[1]+1,\
                            fes_crop_cv3_p[0]:fes_crop_cv3_p[1]+1]\
                                    += gauss[gauss_crop_cv1[0]:gauss_crop_cv1[1]+1,\
                                             gauss_crop_cv2_p[0]:gauss_crop_cv2_p[1]+1,\
                                             gauss_crop_cv3_p[0]:gauss_crop_cv3_p[1]+1]\
                                    * self.heights[line]
                
                if self.periodic[0] and self.periodic[1] and self.periodic[2]:
                    if ((cv1bin[line] < gauss_center_to_end) or (cv1bin[line] > (self.res-gauss_center_to_end)))\
                            and ((cv2bin[line] < gauss_center_to_end) or (cv2bin[line] > (self.res-gauss_center_to_end))) \
                            and ((cv3bin[line] < gauss_center_to_end) or (cv3bin[line] > (self.res-gauss_center_to_end))) :
                        fes[fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1,\
                            fes_crop_cv2_p[0]:fes_crop_cv2_p[1]+1,\
                            fes_crop_cv3_p[0]:fes_crop_cv3_p[1]+1]\
                                    += gauss[gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1,\
                                             gauss_crop_cv2_p[0]:gauss_crop_cv2_p[1]+1,\
                                             gauss_crop_cv3_p[0]:gauss_crop_cv3_p[1]+1]\
                                    * self.heights[line]
            
            if print_output:
                print("\n")
            if subtract_min:
                fes = fes-np.min(fes)
            self.fes = np.array(fes)
        else:
            print("Fes object doesn't have supported number of CVs.")
    
    def makefes2(self, resolution, time_min, time_max, print_output=True, subtract_min = True):
        """
        Function internally used to sum Hills in the same way as Plumed sum_hills. 
        """
        self.res = resolution
        
        if self.cvs == 1:
            fes = np.zeros((self.res))
            progress = 1
            max_progress = self.res ** self.cvs
            
            for x in range(self.res):
                progress += 1
                if print_output and (((progress) % 200 == 0) or (progress == max_progress)):
                    print(f"Constructing free energy surface: {(progress/max_progress):.2%} finished", end="\r")
                
                dist_cv1 = self.cv1[time_min-1:time_max]-(self.cv1min+(x)*self.cv1_fes_range/(self.res))
                if self.periodic[0]:
                    dist_cv1[dist_cv1<-0.5*self.cv1_fes_range] += self.cv1_fes_range
                    dist_cv1[dist_cv1>+0.5*self.cv1_fes_range] -= self.cv1_fes_range
                
                dp2 = dist_cv1**2/(2*self.s1[time_min-1:time_max]**2)
                tmp = np.zeros(time_max-time_min+1)
                heights = self.heights[time_min-1:time_max]
                tmp[dp2<6.25] = heights[dp2<6.25] * (np.exp(-dp2[dp2<6.25]) * 1.00193418799744762399 - 0.00193418799744762399)
                fes[x] = -tmp.sum()
                    
            if subtract_min:
                fes = fes - np.min(fes)
            self.fes = np.array(fes)
            if print_output:
                print("\n")
            
        elif self.cvs == 2:
            fes = np.zeros((self.res, self.res))
            progress = 0
            max_progress = self.res ** self.cvs
            
            for x in range(self.res):
                dist_cv1 = self.cv1[time_min-1:time_max]-(self.cv1min+(x)*self.cv1_fes_range/(self.res))
                if self.periodic[0]:
                    dist_cv1[dist_cv1<-0.5*self.cv1_fes_range] += self.cv1_fes_range
                    dist_cv1[dist_cv1>+0.5*self.cv1_fes_range] -= self.cv1_fes_range
                    
                for y in range(self.res):
                    progress += 1
                    if print_output and ((progress) % 200 == 0 or (progress == max_progress)):
                        print(f"Constructing free energy surface: {(progress/max_progress):.2%} finished", end="\r")
                    
                    dist_cv2 = self.cv2[time_min-1:time_max]-(self.cv2min+(y)*self.cv2_fes_range/(self.res))
                    if self.periodic[1]:
                        dist_cv2[dist_cv2<-0.5*self.cv2_fes_range] += self.cv2_fes_range
                        dist_cv2[dist_cv2>+0.5*self.cv2_fes_range] -= self.cv2_fes_range
                        
                    dp2 = dist_cv1**2/(2*self.s1[time_min-1:time_max]**2) + dist_cv2**2/(2*self.s2[time_min-1:time_max]**2)
                    tmp = np.zeros(time_max-time_min+1)
                    heights = self.heights[time_min-1:time_max]
                    tmp[dp2<6.25] = heights[dp2<6.25] * (np.exp(-dp2[dp2<6.25]) * 1.00193418799744762399 - 0.00193418799744762399)
                    fes[x,y] = -tmp.sum()
                    
            if subtract_min:
                fes = fes - np.min(fes)
            self.fes = np.array(fes)
            if print_output:
                print("\n")
            
        elif self.cvs == 3:
            fes = np.zeros((self.res, self.res, self.res))
            progress = 0
            max_progress = self.res ** self.cvs
            
            for x in range(self.res):
                dist_cv1 = self.cv1[time_min-1:time_max]-(self.cv1min+(x)*self.cv1_fes_range/(self.res))
                if self.periodic[0]:
                    dist_cv1[dist_cv1<-0.5*self.cv1_fes_range] += self.cv1_fes_range
                    dist_cv1[dist_cv1>+0.5*self.cv1_fes_range] -= self.cv1_fes_range
                    
                for y in range(self.res):
                    dist_cv2 = self.cv2[time_min-1:time_max]-(self.cv2min+(y)*self.cv2_fes_range/(self.res))
                    if self.periodic[1]:
                        dist_cv2[dist_cv2<-0.5*self.cv2_fes_range] += self.cv2_fes_range
                        dist_cv2[dist_cv2>+0.5*self.cv2_fes_range] -= self.cv2_fes_range
                        
                    for z in range(self.res):
                        progress += 1
                        if print_output and ((progress) % 200 == 0 or (progress == max_progress)):
                            print(f"Constructing free energy surface: {(progress/max_progress):.2%} finished", end="\r")
                        
                        dist_cv3 = self.cv3[time_min-1:time_max]-(self.cv3min+(z)*self.cv3_fes_range/(self.res))
                        if self.periodic[2]:
                            dist_cv3[dist_cv3<-0.5*self.cv3_fes_range] += self.cv3_fes_range
                            dist_cv3[dist_cv3>+0.5*self.cv3_fes_range] -= self.cv3_fes_range
                        
                        dp2 = dist_cv1**2/(2*self.s1[time_min-1:time_max]**2) + \
                              dist_cv2**2/(2*self.s2[time_min-1:time_max]**2) + \
                              dist_cv3**2/(2*self.s3[time_min-1:time_max]**2)
                        tmp = np.zeros(time_max-time_min+1)
                        heights = self.heights[time_min-1:time_max]
                        tmp[dp2<6.25] = heights[dp2<6.25] * (np.exp(-dp2[dp2<6.25]) * 1.00193418799744762399 - 0.00193418799744762399)
                        fes[x,y,z] = -tmp.sum()
                        
            if subtract_min:
                fes = fes - np.min(fes)
            self.fes = np.array(fes)
            if print_output:
                print("\n")
        else:
            print(f"Error: unsupported number of CVs: {self.cvs}.")
    
    def plot(self, png_name=None, contours=True, contours_spacing=0.0, aspect = 1.0, cmap = "RdYlBu_r", 
             energy_unit="kJ/mol", xlabel=None, ylabel=None, zlabel=None, label_size=12, clabel_size = 12,
             image_size=None, image_size_unit="in", dpi=100, vmin = 0, vmax = None, 
             opacity=0.2, levels=None, title = None, off_screen = False, 
             xlim=[None, None], ylim=[None, None], return_fig=False):
        """
        Function used to visualize FES, based on Matplotlib for 1D and 2D FES or PyVista for 3D FES. 
        
        ```python
        fes.plot(png_name="fes.png")
        ```
        
        Parameters:
        
        * png_name = String. If this parameter is supplied, the picture of FES will be saved under this name to the current working directory. 
        For 1D and 2D FES, the recommended format in '.png'. For 3D FES, the formats supported by PyVista are '.svg','.eps','.ps','.pdf' and '.tex'. 
        
        * contours (default=True) = whether contours should be shown on 2D FES
        
        * contours_spacing (default=0.0) = when a positive number is set, it will be used as spacing for contours on 2D FES. 
                Otherwise, if contours=True, there will be five equally spaced contour levels.
        
        * aspect (default = 1.0) = aspect ratio of the graph. Works with 1D and 2D FES. 
        
        * cmap (default = "RdYlBu_r") = Matplotlib colormap used to color 2D or 3D FES
        
        * energy_unit (default="kJ/mol") = String, used in description of colorbar

        * xlim, ylim (default = [None, None] for both) = list of two values specifying the range of x and y axes respectively. 
        None means that Matplotlib will choose appropriate range, so this keyword is useful when you need to overwrite it. 
        Does not work for 3D FES at the moment. 
        
        * xlabel, ylabel, zlabel = Strings, if provided, they will be used as labels for the graphs
        
        * labelsize (default = 12) = size of text in labels
        
        * image_size (default = [9,6]) = List of the width and height of the picture

        * image_size_unit (default = "in") = Units for width and height of the picture, accepts "in" as inches, "cm", "mm" or "px" as pixels. 

        * dpi (default = 100) = DPI of the resulting image (for 1D and 2D FES). 
        
        * vmin (default=0) = real number, lower bound for the colormap on 2D FES
        
        * vmax = real number, upper bound for the colormap on 2D FES
        
        * opacity (default=0.2) = number between 0 and 1, is the opacity of isosurfaces of 3D FES
        
        * levels = Here you can specify list of free energy values for isosurfaces on 3D FES. 
                        If not provided, default values from contours parameters will be used instead. 

        * title = optional, string that defines the title of the graph

        * offscreen (default = False) = for 3D FES only, grapf will not be shown after creation - used internally by metadynminer when making animations

        * return_fig (default=False) = whether the method should return the Matplotlib.Pyplot.figure object for further use. 
        In the case of plotting 3D FES, it returns Pyvista.plotter object instead. 
        """
        if vmax == None:
            vmax = np.max(self.fes)+0.01 # if the addition is smaller than 0.01, the 3d plot stops working. 
            
        if contours_spacing == 0.0:
            contours_spacing = (vmax-vmin)/5.0
        
        cmap = cm.get_cmap(cmap)
        
        cmap.set_over("white")
        cmap.set_under("white")
        
        if self.cvs >= 1:
            cv1min = self.cv1min
            cv1max = self.cv1max
        if self.cvs >=2:
            cv2min = self.cv2min
            cv2max = self.cv2max 
        if self.cvs == 3:
            cv3min = self.cv3min
            cv3max = self.cv3max
        
        if image_size == None:
            image_size = [9,6]    
        
        if image_size_unit == "cm":
            image_size[0] /= 2.54
            image_size[1] /= 2.54
        elif image_size_unit == "mm":
            image_size[0] /= 25.4
            image_size[1] /= 25.4
        elif image_size_unit == "px":
            image_size[0] /= dpi
            image_size[1] /= dpi
        elif image_size_unit != "in":
            print(f"Warning: unknown image_size_unit value: {image_size_unit}. Using inches instead. ")
            
        
        if self.cvs == 1:
            fig = plt.figure(figsize=(image_size[0],image_size[1]), dpi=dpi)
            X = np.linspace(cv1min, cv1max, self.res)
            plt.plot(X, self.fes)
            if xlabel == None:
                plt.xlabel(f'CV1 - {self.cv1_name}', size=label_size)
            else:
                plt.xlabel(xlabel, size=label_size)
            if ylabel == None:
                plt.ylabel(f'free energy ({energy_unit})', size=label_size)
            else:
                plt.ylabel(ylabel, size=label_size)
        
            if title != None:
                plt.title(title, size=label_size)
    
            plt.xticks(fontsize=label_size)
            plt.yticks(fontsize=label_size)
                
            ax = plt.gca()
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            
            if png_name != None:
                plt.savefig(png_name, bbox_inches = 'tight')

            if return_fig:
                return fig
            
        if self.cvs == 2:
            fig = plt.figure(figsize=(image_size[0], image_size[1]), dpi=dpi)
            plt.imshow(np.rot90(self.fes, axes=(0,1)), cmap=cmap, interpolation='nearest', 
                       extent=[cv1min, cv1max, cv2min, cv2max], 
                       aspect = (((cv1max-cv1min)/(cv2max-cv2min))/(aspect)),
                       vmin = vmin, vmax = vmax)
            cbar = plt.colorbar()
            cbar.ax.tick_params(labelsize=label_size) 
            cbar.set_label(energy_unit, size=label_size)
            
            if contours:
                if levels != None:
                    cont = plt.contour(np.rot90(self.fes, axes=(0,1)), 
                             levels = levels, 
                             extent=[cv1min, cv1max, cv2max, cv2min], 
                             colors = "k", linestyles = "solid")
                    plt.clabel(cont, levels = levels, fontsize=clabel_size)
                else:
                    cont = plt.contour(np.rot90(self.fes, axes=(0,1)), 
                             levels = np.arange(vmin, (vmax - 0.01), contours_spacing), 
                             extent=[cv1min, cv1max, cv2max, cv2min], 
                             colors = "k", linestyles = "solid")
                    plt.clabel(cont, levels = np.arange(vmin, (vmax - 0.01), contours_spacing), fontsize=clabel_size)
            
            if xlabel == None:
                plt.xlabel(f'CV1 - {self.cv1_name}', size=label_size)
            else:
                plt.xlabel(xlabel, size=label_size)
                
            if ylabel == None:
                plt.ylabel(f'CV2 - {self.cv2_name}', size=label_size)
            else:
                plt.ylabel(ylabel, size=label_size)
            
            if title != None:
                plt.title(title, size=label_size)
    
            plt.xticks(fontsize=label_size)
            plt.yticks(fontsize=label_size)

            ax = plt.gca()
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            ax.tick_params(axis='x', labelsize=label_size)
            ax.tick_params(axis='y', labelsize=label_size)
            
            if png_name != None:
                plt.savefig(png_name, bbox_inches = 'tight')
            
            if return_fig:
                return fig
        
        if self.cvs == 3:
            if xlabel == None:
                xlabel = "CV1 - " + self.cv1_name
            if ylabel == None:
                ylabel = "CV2 - " + self.cv2_name
            if zlabel == None:
                zlabel = "CV3 - " + self.cv3_name
            
            grid = pv.ImageData(
                dimensions=(self.res, self.res, self.res),
                spacing=((cv1max-cv1min)/self.res,(cv2max-cv2min)/self.res,(cv3max-cv3min)/self.res),
                origin=(cv1min, cv2min, cv3min)
            )
            grid["vol"] = self.fes.ravel(order="F")
            if levels == None:
                contours = grid.contour(np.arange(0, (vmax - 0.01), contours_spacing))
            else:
                contours = grid.contour(levels)
            fescolors = []
            for i in range(contours.points.shape[0]):
                fescolors.append(self.fes[int((contours.points[i,0]-cv1min)*self.res/(cv1max-cv1min)),
                                          int((contours.points[i,1]-cv2min)*self.res/(cv2max-cv2min)),
                                          int((contours.points[i,2]-cv3min)*self.res/(cv3max-cv3min))])
            #%% Visualization
            bounds = [cv1min, cv1max, cv2min, cv2max, cv3min, cv3max]
            pv.set_plot_theme('document')
            p = pv.Plotter()
            p.add_mesh(contours, scalars=fescolors, opacity=opacity, cmap=cmap, show_scalar_bar=False, interpolate_before_map=True)
            p.show_bounds(bounds=bounds, xtitle=xlabel, ytitle=ylabel, ztitle=zlabel, grid=True, location="outer")
            if title != None:
                p.add_title(title, font_size=label_size)
            if not off_screen:
                p.show()
            
            if png_name != None:
                p.save_graphic(png_name)
            
            if return_fig:
                return p
    
    def set_fes(self, fes):
        self.fes = fes
        
    def surface_plot(self, png_name=None, cmap = "RdYlBu_r", 
                     energy_unit="kJ/mol", xlabel=None, ylabel=None, zlabel=None, dpi=100, label_size=12, image_size=None, image_size_unit="in", rstride=1, cstride=1, vmin = 0, vmax = None, title=None, return_fig=False):
        """
        Function for visualization of 2D FES as 3D surface plot. For now, it is based on Matplotlib, but there are issues with interactivity. 
        
        It can be interacted with in jupyter notebook or jupyter lab in %matplotlib widget mode. Otherwise it is just static image of the 3D surface plot. 
        
        ```python
        %matplotlib widget
        fes.surface_plot()
        ```
        
        """
        if self.cvs == 2:
            cv1min = self.cv1min
            cv1max = self.cv1max
            cv2min = self.cv2min
            cv2max = self.cv2max 
            
            x = np.linspace(cv1min, cv1max, self.res)
            y = np.linspace(cv2min, cv2max, self.res)
            
            X, Y = np.meshgrid(x, y)
            Z = self.fes.T
    
            if image_size == None:
                image_size = [9,9]
            
            if image_size_unit == "cm":
                image_size[0] /= 2.54
                image_size[1] /= 2.54
            elif image_size_unit == "mm":
                image_size[0] /= 25.4
                image_size[1] /= 25.4
            elif image_size_unit == "px":
                image_size[0] /= dpi
                image_size[1] /= dpi
            elif image_size_unit != "in":
                print(f"Warning: unknown image_size_unit value: {image_size_unit}. Using inches instead. ")
            
            fig = plt.figure(figsize=(image_size[0],image_size[1]), dpi=dpi)
            ax = fig.add_subplot(projection="3d")
            ax.plot_surface(X,Y,Z, cmap=cmap, rstride=rstride, cstride=cstride)
            
            if xlabel == None:
                ax.set_xlabel(f'CV1 - {self.cv1_name}', size=label_size)
            else:
                ax.set_xlabel(xlabel, size=label_size)
            if ylabel == None:
                ax.set_ylabel(f'CV2 - {self.cv2_name}', size=label_size)
            else:
                ax.set_ylabel(ylabel, size=label_size)
            if zlabel == None:
                ax.set_zlabel(f'free energy ({energy_unit})', size=label_size)
            else:
                ax.set_zlabel(zlabel, size=label_size)
    
            if title != None:
                plt.title(title, size=label_size)
    
            plt.xticks(fontsize=label_size)
            plt.yticks(fontsize=label_size)
    
            if png_name != None:
                plt.savefig(png_name, bbox_inches = 'tight')
            
            if return_fig:
                return fig
                
        else:
            print(f"Error: Surface plot only works for FES with exactly two CVs, and this FES has {self.cvs} CV.")
    
    def remove_CV(self, CV=None, energy_unit="kJ/mol", temp=300.0):
        """
        This function is used to remove a CV from an existing FES. The function first recalculates the FES to an array of probabilities. The probabilities 
        are summed along the CV to be removed, and resulting probability distribution with 1 less dimension is converted back to FES. 

        ```python
        fes_cv1 = fes.remove_CV(CV=2)
        ```
        
        Parameters:
        
        * CV = integer, the CV to be removed
        
        * energy_unit (default="kJ/mol") = has to be either "kJ/mol" or "kcal/mol". Make sure to suply the correct energy unit, otherwise you will get wrong FES as a result. 
        
        * temp (default=300.0) = temperature of the simulation in Kelvins.
        """
        CV = int(float(CV))
        print(f"Removing CV {CV}.")
        if CV > self.cvs:
            print("Error: The CV to remove is not available in this FES object.")
            return None
        if self.cvs == 1:
            print("Error: You can not remove the only CV. ")
            return None
        elif self.cvs == 2:
            if energy_unit == "kJ/mol":
                probabilities = np.exp(-1000*self.fes/8.314/temp)
                if CV == 1:
                    new_prob = np.sum(probabilities, axis=0)
                    new_fes = Fes(hills=None)
                    new_fes.fes = -8.314*temp*np.log(new_prob)/1000
                    new_fes.fes = new_fes.fes - np.min(new_fes.fes)
                    new_fes.cvs = 1
                    new_fes.res = self.res
                    new_fes.periodic = [self.periodic[1]]
                    new_fes.cv1min = self.cv2min
                    new_fes.cv1max = self.cv2max
                    new_fes.cv1_name = self.cv2_name
                    new_fes.cv1per = self.cv2per
                    new_fes.cv1_fes_range = self.cv2_fes_range
                if CV == 2:
                    new_prob = np.sum(probabilities, axis=1)
                    new_fes = Fes(hills=None)
                    new_fes.fes = -8.314*temp*np.log(new_prob)/1000
                    new_fes.fes = new_fes.fes - np.min(new_fes.fes)
                    new_fes.cvs = 1
                    new_fes.res = self.res
                    new_fes.periodic = [self.periodic[0]]
                    new_fes.cv1min = self.cv1min
                    new_fes.cv1max = self.cv1max
                    new_fes.cv1_name = self.cv1_name
                    new_fes.cv1per = self.cv1per
                    new_fes.cv1_fes_range = self.cv1_fes_range
                return new_fes
            elif energy_unit == "kcal/mol":
                probabilities = np.exp(-1000*4.184*self.fes/8.314/temp)
                if CV == 1:
                    new_prob = np.sum(probabilities, axis=1)
                    new_fes = Fes(hills=None)
                    new_fes.fes = -8.314*temp*np.log(new_prob)/1000/4.184
                    new_fes.fes = new_fes.fes - np.min(new_fes.fes)
                    new_fes.cvs = 1
                    new_fes.res = self.res
                    new_fes.periodic = [self.periodic[1]]
                    new_fes.cv1min = self.cv2min
                    new_fes.cv1max = self.cv2max
                    new_fes.cv1_name = self.cv2_name
                    new_fes.cv1per = self.cv2per
                    new_fes.cv1_fes_range = self.cv2_fes_range
                if CV == 2:
                    new_prob = np.sum(probabilities, axis=0)
                    new_fes = Fes(hills=None)
                    new_fes.fes = -8.314*temp*np.log(new_prob)/1000/4.184
                    new_fes.fes = new_fes.fes - np.min(new_fes.fes)
                    new_fes.cvs = 1
                    new_fes.res = self.res
                    new_fes.periodic = [self.periodic[0]]
                    new_fes.cv1min = self.cv1min
                    new_fes.cv1max = self.cv1max
                    new_fes.cv1_name = self.cv1_name
                    new_fes.cv1per = self.cv1per
                    new_fes.cv1_fes_range = self.cv1_fes_range
                return new_fes
            else:
                print("Error: unknown energy unit")
                return None
        elif self.cvs == 3:
            if energy_unit == "kJ/mol":
                probabilities = np.exp(-1000*self.fes/8.314/temp)
                if CV == 1:
                    new_prob = np.sum(probabilities, axis=0)
                    new_fes = Fes(hills=None)
                    new_fes.fes = -8.314*temp*np.log(new_prob)/1000
                    new_fes.fes = new_fes.fes - np.min(new_fes.fes)
                    new_fes.cvs = 2
                    new_fes.res = self.res
                    new_fes.periodic = [self.periodic[1], self.periodic[2]]
                    new_fes.cv1min = self.cv2min
                    new_fes.cv1max = self.cv2max
                    new_fes.cv2min = self.cv3min
                    new_fes.cv2max = self.cv3max
                    new_fes.cv1_name = self.cv2_name
                    new_fes.cv2_name = self.cv3_name
                    new_fes.cv1per = self.cv2per
                    new_fes.cv2per = self.cv3per
                    new_fes.cv1_fes_range = self.cv2_fes_range
                    new_fes.cv2_fes_range = self.cv3_fes_range
                if CV == 2:
                    new_prob = np.sum(probabilities, axis=1)
                    new_fes = Fes(hills=None)
                    new_fes.fes = -8.314*temp*np.log(new_prob)/1000
                    new_fes.fes = new_fes.fes - np.min(new_fes.fes)
                    new_fes.cvs = 2
                    new_fes.res = self.res
                    new_fes.periodic = [self.periodic[0], self.periodic[2]]
                    new_fes.cv1min = self.cv1min
                    new_fes.cv1max = self.cv1max
                    new_fes.cv2min = self.cv3min
                    new_fes.cv2max = self.cv3max
                    new_fes.cv1_name = self.cv1_name
                    new_fes.cv2_name = self.cv3_name
                    new_fes.cv1per = self.cv1per
                    new_fes.cv2per = self.cv3per
                    new_fes.cv1_fes_range = self.cv1_fes_range
                    new_fes.cv2_fes_range = self.cv3_fes_range
                if CV == 3:
                    new_prob = np.sum(probabilities, axis=2)
                    new_fes = Fes(hills=None)
                    new_fes.fes = -8.314*temp*np.log(new_prob)/1000
                    new_fes.fes = new_fes.fes - np.min(new_fes.fes)
                    new_fes.cvs = 2
                    new_fes.res = self.res
                    new_fes.periodic = [self.periodic[0], self.periodic[1]]
                    new_fes.cv1min = self.cv1min
                    new_fes.cv1max = self.cv1max
                    new_fes.cv2min = self.cv2min
                    new_fes.cv2max = self.cv2max
                    new_fes.cv1_name = self.cv1_name
                    new_fes.cv2_name = self.cv2_name
                    new_fes.cv1per = self.cv1per
                    new_fes.cv2per = self.cv2per
                    new_fes.cv1_fes_range = self.cv1_fes_range
                    new_fes.cv2_fes_range = self.cv2_fes_range
                return new_fes
            elif energy_unit == "kcal/mol":
                probabilities = np.exp(-1000*4.184*self.fes/8.314/temp)
                if CV == 1:
                    new_prob = np.sum(probabilities, axis=0)
                    new_fes = Fes(hills=None)
                    new_fes.fes = -8.314*temp*np.log(new_prob)/1000/4.184
                    new_fes.fes = new_fes.fes - np.min(new_fes.fes)
                    new_fes.cvs = 2
                    new_fes.res = self.res
                    new_fes.periodic = [self.periodic[1], self.periodic[2]]
                    new_fes.cv1min = self.cv2min
                    new_fes.cv1max = self.cv2max
                    new_fes.cv2min = self.cv3min
                    new_fes.cv2max = self.cv3max
                    new_fes.cv1_name = self.cv2_name
                    new_fes.cv2_name = self.cv3_name
                    new_fes.cv1per = self.cv2per
                    new_fes.cv2per = self.cv3per
                    new_fes.cv1_fes_range = self.cv2_fes_range
                    new_fes.cv2_fes_range = self.cv3_fes_range
                if CV == 2:
                    new_prob = np.sum(probabilities, axis=1)
                    new_fes = Fes(hills=None)
                    new_fes.fes = -8.314*temp*np.log(new_prob)/1000/4.184
                    new_fes.fes = new_fes.fes - np.min(new_fes.fes)
                    new_fes.cvs = 2
                    new_fes.res = self.res
                    new_fes.periodic = [self.periodic[0], self.periodic[2]]
                    new_fes.cv1min = self.cv1min
                    new_fes.cv1max = self.cv1max
                    new_fes.cv2min = self.cv3min
                    new_fes.cv2max = self.cv3max
                    new_fes.cv1_name = self.cv1_name
                    new_fes.cv2_name = self.cv3_name
                    new_fes.cv1per = self.cv1per
                    new_fes.cv2per = self.cv3per
                    new_fes.cv1_fes_range = self.cv1_fes_range
                    new_fes.cv2_fes_range = self.cv3_fes_range
                if CV == 3:
                    new_prob = np.sum(probabilities, axis=2)
                    new_fes = Fes(hills=None)
                    new_fes.fes = -8.314*temp*np.log(new_prob)/1000/4.184
                    new_fes.fes = new_fes.fes - np.min(new_fes.fes)
                    new_fes.cvs = 2
                    new_fes.res = self.res
                    new_fes.periodic = [self.periodic[0], self.periodic[1]]
                    new_fes.cv1min = self.cv1min
                    new_fes.cv1max = self.cv1max
                    new_fes.cv2min = self.cv2min
                    new_fes.cv2max = self.cv2max
                    new_fes.cv1_name = self.cv1_name
                    new_fes.cv2_name = self.cv2_name
                    new_fes.cv1per = self.cv1per
                    new_fes.cv2per = self.cv2per
                    new_fes.cv1_fes_range = self.cv1_fes_range
                    new_fes.cv2_fes_range = self.cv2_fes_range
                return new_fes
            else:
                print("Error: unknown energy unit")
                return None

    def flooding_animation(self, gif_name = "flooding.gif", use_vmax_from_end = True, with_minima = True, use_minima_from_end=False, cmap="RdYlBu_r", xlabel=None, ylabel=None, zlabel=None, label_size=12, image_size=None, image_size_unit="in", dpi=100, tu = "ps", time_min=None, time_max=None, step=1000, contours_spacing = 20, levels=None, opacity = 0.2, vmin = 0, vmax = None, energy_unit="kJ/mol", clear_temporary_folder=True, temporary_folder_name="temporary_folder", time_unit="ps", fps=5, enable_loop = True, minima_precise = False, minima_nbins=8, temp=300.0, 
                           xlim=[None, None], ylim=[None, None]):
        """
        This method is used to make an animation that shows, how the FES was evolving during metadynamics simulation. It creates temporary folder and svaes plots of FES at different times during simulation, then it concatenates them to make a gif animation and removes the temporary files (remove can be switched off, if necessary). 

        Use:
        
        ```python
        fes.flooding_animation()
        ```
        Parameters:
        * gif_name (default="flooding.gif") = name for the animation file

        * optional parameters used when calling
        ```python
        fes.plot()
        ```
        or
        ```python
        minima.plot()
        ```
        can be used to adjust the graphs used to create animation

        * use_vmax_from_end (default=True), boolean value, if it is True and vmax is not specified, it will use the vmax value from the end point of the fes object for all frames of animation; then, for 2D FES, the colors will be comparable to each other

        * with_minima (default = True), boolean value, if True, graphs will be shown with letters at each local minima found at each frame of the animation

        * use _minima_from_end (default = True), boolean value, if True, the local minima from the end point of the fes object wil be depicted at each frame of animation

        * minima_precise (default=False) = whether the local minima should be calculated using the precise algorithm (see Minima for more info)

        * minima_nbins (default = 8) = the nbins keyword for minima localization

        * temp (default = 300.0), energy_unit (default = "kJ/mol") - keywords for minima localization in case minima_precise = True

        * step (default=1000), integer, frames for animation will be made at each n-th line in HILLS file

        * temporary_folder_name (default="temporary_folder"), name of the temporary folder where the individual graphs will be saved; directory with this name shouldn't be present in working irectory when calling the method, otherwise it will throw an error
        
        * clear_temporary_folder (default=True), if set to false, the graphs made for each frame of the simulation will not be removed afterwards, so they can be viewed later

        * fps (defalt=5) = how many frames per second the animation will have

        * enable_loop (default=True) = whether the animation should be running in loop
        """
        tu = TU(tu)
        
        current_directory = os.getcwd()
        final_directory = os.path.join(current_directory, temporary_folder_name)
        if os.path.exists(final_directory):
            print(f"Error: directory with the name {temporary_folder_name} already exists. Try using different name using temporary_folder_name keyword. ")
            return None
        else:
           os.makedirs(final_directory)
        flooding_fes = copy.deepcopy(self)
        step_fes = copy.deepcopy(self)
        
        if time_min == time_max == 0:
            print("Error: Values of start and end time are zero.")
            return None
        if time_min != None:
            time_min = tu.inps(time_min)
            if time_min < 0:
                print("Warning: Start time is lower than zero, it will be set to zero instead. ")
                time_min = 0
            if time_min < int(self.hills.hills[0,0]):
                print(f"Warning: Start time {tu.inps(time_min)} ps is lower than the first time from HILLS file {int(self.hills.hills[0,0])}, which will be used instead. ")
                time_min = int(self.hills.hills[0,0])
        else:
            time_min = int(self.hills.hills[0,0])
        if time_max != None:
            time_max = tu.inps(time_max)
            if time_max < time_min:
                print("Warning: End time is lower than start time. Values are flipped. ")
                time_value = time_max
                time_max = time_min
                time_min = time_value
            if time_max > int(self.hills.hills[-1,0]):
                print(f"Warning: End time {tu.inps(time_max)} ps is higher than number of lines in HILLS file {int(self.hills.hills[-1,0])}, which will be used instead. ")
                time_max = int(self.hills.hills[-1,0])
        else:
            time_max = int(self.hills.hills[-1,0])
        
        if not self.ignoretime:
            time_max = int(round(((time_max - self.hills.hills[0,0])/self.dt),0)) + 1
            time_min = int(round(((time_min - self.hills.hills[0,0])/self.dt),0)) + 1

        if time_min==None:
            time_min=1
        if time_max==None:
            time_max = int(self.hills.hills[-1,0])

        if (vmax == None) and use_vmax_from_end:
            vmax = np.max(self.fes)+0.1
        if self.cvs == 1 or self.cvs == 2:
            suffix="png"
        if self.cvs == 3:
            suffix="eps"
        
        times = np.array((range(time_min-1, time_max+1, step)))
        image_files = [f'{final_directory}/{times[i]}.{suffix}'.format(i) for i in range(1, len(times))]
        
        minima_final = Minima(self, precise=minima_precise, nbins=minima_nbins, energy_unit=energy_unit, temp=temp, print_output=False)
        
        for i in range(1,len(times)):
            print(f"Constructing flooding animation: {((i+1)/len(times)):.2%} finished", end="\r")
            if self.original==False:
                step_fes.makefes(flooding_fes.res, time_min=(times[i-1]+1), time_max=times[i], print_output=False)
            else:
                step_fes.makefes2(flooding_fes.res, time_min=(times[i-1]+1), time_max=times[i], print_output=False)
            if i == 1:
                flooding_fes.fes = step_fes.fes
            else:
                flooding_fes.fes += step_fes.fes
            flooding_fes.fes = flooding_fes.fes - np.min(flooding_fes.fes)
            try:
                if with_minima:
                    mf = Minima(flooding_fes, precise=minima_precise, print_output=False, nbins=minima_nbins, temp=temp, energy_unit=energy_unit)
                    if use_minima_from_end:
                        mf.minima = minima_final.minima
                    mf.plot(contours_spacing=contours_spacing, cmap = cmap, 
                            xlabel = xlabel, ylabel = ylabel, zlabel=zlabel, 
                            label_size=label_size, image_size=image_size, 
                            image_size_unit=image_size_unit, dpi=dpi, vmin = vmin, 
                            vmax = vmax, levels=levels, energy_unit=energy_unit, 
                            off_screen = True, opacity=opacity, 
                            png_name=f"{final_directory}/{times[i]}.{suffix}", title=f"{times[i]} {time_unit}", xlim=xlim, ylim=ylim)
                else:
                    flooding_fes.plot(contours_spacing=contours_spacing, cmap = cmap, xlabel = xlabel, ylabel = ylabel, zlabel=zlabel, label_size=label_size, image_size=image_size, image_size_unit=image_size_unit, dpi = dpi, vmin = vmin, vmax = vmax, levels=levels, energy_unit=energy_unit, off_screen = True, opacity=opacity, png_name=f"{final_directory}/{times[i]}.{suffix}", title=f"{times[i]} {time_unit}", xlim=xlim, ylim=ylim)
                plt.close()
            except ValueError:
                print("Warning: The first frame of animation would be blank with the current settings, but PyVista 3D plotter can not plot empty meshes. Try to increase the timestep between frames or decrease the spacing between isosurfaces.")
        print("\n")

        duration = 1000/fps
        if enable_loop:
            with imageio.get_writer(gif_name, format="GIF", duration=duration, loop=0) as writer:
                for image_file in image_files:
                    try:
                        image = imageio.imread(image_file)
                        writer.append_data(image)
                    except FileNotFoundError:
                        print("Warning: File for animation was not found.")
        else:
            with imageio.get_writer(gif_name, format="GIF", fps = fps) as writer:
                for image_file in image_files:
                    try:
                        image = imageio.imread(image_file)
                        writer.append_data(image)
                    except FileNotFoundError:
                        print("Warning: File for animation was not found.")

        if clear_temporary_folder:
            for image_file in image_files:
                try:
                    os.remove(image_file)
                except FileNotFoundError:
                        print("Warning: File for animation was not found.")
            try: 
                os.rmdir(final_directory)
            except OSError:
                print(f"Warning: directory{final_directory} is not empty and will not be removed. Metadynminer's temporary files inside were removed. ")
               

    
    def make_gif(self, gif_name=None, cmap = "RdYlBu_r", energy_unit="kJ/mol",
                 xlabel=None, ylabel=None, zlabel=None, label_size=12, image_size=[10,7], 
                  opacity=0.2, levels=None, frames=64):
        """
        Function that generates animation of 3D FES showing different isosurfaces.
        
        ```python
        fes.make_gif(gif_name="FES_animation.gif")
        ```
        
        Parameters:
        
        * gif_name (default="FES.gif") = String. Name of the gif of FES that will be saved in the current working directory.
        
        * cmap (default = "RdYlBu_r") = Matplotlib colormap used to color the 3D FES
        
        * xlabel, ylabel, zlabel = Strings, if provided, they will be used as labels for the graph
        
        * labelsize (default = 12) = size of text in labels
        
        * image_size (default = [10,7]) = List of the width and height of the picture
        
        * opacity (default = 0.2) = number between 0 and 1, is the opacity of isosurfaces of 3D FES
        
        * levels = Here you can specify list of free energy values for isosurfaces on 3D FES. 
                If not provided, default values from contours parameters will be used instead. 
        
        * frames (default = 64) = Number of frames the animation will be made of. 
        """
        if self.cvs == 3:
            cv1min = self.cv1min
            cv1max = self.cv1max
            cv2min = self.cv2min
            cv2max = self.cv2max
            cv3min = self.cv3min
            cv3max = self.cv3max
            
            values = np.linspace(np.min(self.fes)+0.1, np.max(self.fes), num=frames)
            grid = pv.ImageData(
                dimensions=(self.res, self.res, self.res),
                spacing=((cv1max-cv1min)/self.res,(cv2max-cv2min)/self.res,(cv3max-cv3min)/self.res),
                origin=(cv1min, cv2min, cv3min),
            )
            grid["vol"] = self.fes.ravel(order="F")
            surface = grid.contour(values[:1])
            surfaces = [grid.contour([v]) for v in values]
            surface = surfaces[0].copy()
            
            pv.set_plot_theme('document')
            plotter = pv.Plotter(off_screen=True)
            # Open a movie file
            if gif_name == None:
                gif_name = "FES_animation.gif"
            plotter.open_gif(gif_name)

            # Add initial mesh
            plotter.add_mesh(
                surface,
                opacity=opacity,
                clim=grid.get_data_range(),
                show_scalar_bar=False,
                cmap="RdYlBu_r"
            )
            plotter.add_mesh(grid.outline_corners(), color="k")
            text = plotter.add_text(f"{values[0]:.2f}+kJ/mol", position='lower_right', font_size=12)
            if xlabel == None and ylabel == None and zlabel == None:
                plotter.show_grid(xtitle=f"CV1 - {self.cv1_name}", ytitle=f"CV2 - {self.cv2_name}", ztitle=f"CV3 - {self.cv3_name}")
            else:
                plotter.show_grid(xtitle=xlabel, ytitle=ylabel, ztitle=zlabel)
            plotter.set_background('white')
            plotter.show(auto_close=False)

            # Run through each frame
            for surf in range(len(surfaces)):
                surface.copy_from(surfaces[surf])
                plotter.remove_actor(text)
                text = plotter.add_text(f"{values[surf]:.2f} {energy_unit}", position='lower_right', font_size=12)
                plotter.write_frame()  # Write this frame
            # Run through backwards
            for surf in range(len(surfaces)-1,0,-1):
                surface.copy_from(surfaces[surf])
                plotter.remove_actor(text)
                text = plotter.add_text(f"{values[surf]:.2f} {energy_unit}", position='lower_right', font_size=12)
                plotter.write_frame()  # Write this frame

            # Be sure to close the plotter when finished
            plotter.close()
        else:
            print("Error: this method is only available for FES with 3 CVs.")     
 
class Minima():
    """
    Object of Minima class is created to find local free energy minima on FES and calculates their free energy values. 
    
    The list of minima is stored as pandas dataframe. 
    
    Command:
    ```python
    minima = metadynminer.Minima(fes=f, nbins=8)
    ```
    
    List of minima can be later shown like this:
    
    ```python
    print(minima.minima)
    ```
    
    Parameters:
    
    * fes = Fes object to find the minima on
    
    * nbins (default = 8) = number of bins to divide the FES

    * precise (default=True) = if True, the local minima will use an algorithm which finds all local minima, even very shallow and probably unimportant minima, each point on the FES will be assigned to the minimum the system would most likely go to, if it only follows the gradient of free energy, and free energy value of minima will be calculated from each point on FES assigned to the respective minima. This results in more precise free energy values, as it accounts for the width of the minimum as well. For this calculation the unit of free energy and the thermodynamical temperature of the simulation must be supplied. This algorithm doesn't use the nbins keyword. 
    If you set precise = False, the method will use the original algorithm from the metadynminer package for R. In this algorithm the FES is 
    first divided to number of bins (can be set with option nbins, default is 8)
    and the absolute minima is found for each bin. Then the algorithm checks 
    if this point is really a local minimum by comparing to the surrounding points of FES.
    This algorithm only accounts for the depth of each minima, which is less precise, but usually sufficient. 
    In some cases, this algorithm is the prefered one, 
    because on some free energy landscapes the total number of local free energy minima can reach tens of thousands, 
    which makes the calculation using precise algorithm slow and impractical. 

    * temp (default = 300.0) = thermodynamical temperature of the simulation. Used only if precise=True. 

    * energy_unit (default = "kJ/mol") = energy unit of the free energy surface; must be either "kJ/mol" or "kcal/mol". Used only if precise=True. 

    * max_iteration (default=10000), the maximum number of iteration the algorithm will last when assigning FES points to their respective local minima

    Attributes: 

    Minima.fes = source free energy surface 

    Minima.minima = pandas dataframe with local minima and their properties

    Minima.m_fes = if you use precise=True, this attribute contains the FES, but each bin contains the index of the minima to which it belongs

    * and other attributes analogous to those of Fes objects
    """
    
    def __init__(self, fes, nbins = 8, precise=True, temp=300.0, energy_unit="kJ/mol", max_iteration=10000, print_output=True):
        self.fes = fes.fes
        self.periodic = fes.periodic
        self.cvs = fes.cvs
        self.res = fes.res

        if self.cvs >= 1:
            self.cv1_name = fes.cv1_name
            self.cv1min = fes.cv1min
            self.cv1max = fes.cv1max
            self.cv1per = fes.cv1per
            self.cv1_fes_range = fes.cv1_fes_range
        if self.cvs >= 2:
            self.cv2min = fes.cv2min
            self.cv2max = fes.cv2max
            self.cv2_name = fes.cv2_name
            self.cv2per = fes.cv2per
            self.cv2_fes_range = fes.cv2_fes_range
        if self.cvs == 3:
            self.cv3min = fes.cv3min
            self.cv3max = fes.cv3max
            self.cv3_name = fes.cv3_name
            self.cv3per = fes.cv3per
            self.cv3_fes_range = fes.cv3_fes_range

        if precise:
            if energy_unit in ["kJ/mol", "kcal/mol"]:
                self.findminima2(temp=temp, energy_unit=energy_unit, max_iteration=max_iteration, print_output=print_output)
            else:
                print("Error: energy_unit must be either 'kJ/mol' or 'kcal/mol'. ")
        else:
            self.findminima(nbins=nbins)

    def _get_indexes(self, i=0, j=0, k=0):
        # returns array of indexes of the surroundings points of FES, respecting periodicity
        periodic = self.periodic
        if self.cvs == 1:
            a_indexes = np.empty((3))
        
            a_indexes[0] = np.array((i-1)) 
            a_indexes[1] = np.array((i))
            a_indexes[2] = np.array((i+1))
            
            if periodic[0]:
                if i == 0:
                    a_indexes[0] = (self.fes.shape[0]-1)
                elif i == self.fes.shape[0]-1:
                    a_indexes[2] = (0)
            
            # remove uninitialized parts in a_indexes - edges not periodic
            if i == 0 and not periodic[0]:
                a_indexes = a_indexes[1:]
            if i == self.fes.shape[0]-1 and not periodic[0]:
                a_indexes = a_indexes[:-1]
        
            return a_indexes.astype(int)
            
        elif self.cvs == 2:
            a_indexes = np.empty((3,3,2))
        
            a_indexes[0,:,0] = np.array((i-1,i-1,i-1)) 
            a_indexes[1,:,0] = np.array((i,i,i))
            a_indexes[2,:,0] = np.array((i+1,i+1,i+1))
            if i == 0:
                if periodic[0]:
                    a_indexes[0,:,0] = np.array((self.fes.shape[0]-1,self.fes.shape[0]-1,self.fes.shape[0]-1))
            elif i == self.fes.shape[0]-1:
                if periodic[0]:
                    a_indexes[2,:,0] = np.array((0,0,0))
        
            a_indexes[:,0,1] = np.array((j-1,j-1,j-1))
            a_indexes[:,1,1] = np.array((j,j,j))
            a_indexes[:,2,1] = np.array((j+1,j+1,j+1))
            
            if j == 0:
                if periodic[1]:
                    a_indexes[:,0,1] = np.array((self.fes.shape[1]-1,self.fes.shape[1]-1,self.fes.shape[1]-1))
            elif j == self.fes.shape[1]-1:
                if periodic[1]:
                    a_indexes[:,2,1] = np.array((0,0,0))
            
            # remove uninitialized parts in a_indexes - edges not periodic
            if i == 0 and not periodic[0]:
                a_indexes = a_indexes[1:,:,:]
            if i == self.fes.shape[0]-1 and not periodic[0]:
                a_indexes = a_indexes[:-1,:,:]
            if j == 0 and not periodic[1]:
                a_indexes = a_indexes[:,1:,:]
            if j == self.fes.shape[0]-1 and not periodic[1]:
                a_indexes = a_indexes[:,:-1,:]
        
            return a_indexes.astype(int)
        
        elif self.cvs == 3:
            a_indexes = np.empty((3,3,3,3))
        
            a_indexes[0,:,:,0] = np.ones((3,3)) * i-1
            a_indexes[1,:,:,0] = np.ones((3,3)) * i
            a_indexes[2,:,:,0] = np.ones((3,3)) * i+1
            if periodic[0]:
                if i == 0:
                    a_indexes[0,:,:,0] = np.ones((3,3)) * (self.fes.shape[0]-1)
                elif i == self.fes.shape[0]-1:
                    a_indexes[2,:,:,0] = np.zeros((3,3))
        
            a_indexes[:,0,:,1] = np.ones((3,3)) * j-1
            a_indexes[:,1,:,1] = np.ones((3,3)) * j
            a_indexes[:,2,:,1] = np.ones((3,3)) * j+1
            
            if j == 0:
                if periodic[1]:
                    a_indexes[:,0,:,1] = np.ones((3,3)) * (self.fes.shape[1]-1)
            elif j == self.fes.shape[1]-1:
                if periodic[1]:
                    a_indexes[:,2,:,1] = np.zeros((3,3))

            a_indexes[:,:,0,2] = np.ones((3,3)) * k-1
            a_indexes[:,:,1,2] = np.ones((3,3)) * k
            a_indexes[:,:,2,2] = np.ones((3,3)) * k+1
            
            if k == 0:
                if periodic[2]:
                    a_indexes[:,:,0,2] = np.ones((3,3)) * (self.fes.shape[2]-1)
            elif k == self.fes.shape[2]-1:
                if periodic[2]:
                    a_indexes[:,:,2,2] = np.zeros((3,3))
            
            # remove uninitialized parts in a_indexes - edges not periodic
            if i == 0 and not periodic[0]:
                a_indexes = a_indexes[1:,:,:,:]
            if i == self.fes.shape[0]-1 and not periodic[0]:
                a_indexes = a_indexes[:-1,:,:,:]
            if j == 0 and not periodic[1]:
                a_indexes = a_indexes[:,1:,:,:]
            if j == self.fes.shape[1]-1 and not periodic[1]:
                a_indexes = a_indexes[:,:-1,:,:]
            if k == 0 and not periodic[2]:
                a_indexes = a_indexes[:,:,1:,:]
            if k == self.fes.shape[2]-1 and not periodic[2]:
                a_indexes = a_indexes[:,:,:-1,:]
        return a_indexes
    
    def findminima2(self, temp=300.0, energy_unit="kJ/mol", max_iteration=10000, print_output=True):
        if self.cvs >= 1:
            cv1min = self.cv1min
            cv1max = self.cv1max
        if self.cvs >=2:
            cv2min = self.cv2min
            cv2max = self.cv2max 
        if self.cvs == 3:
            cv3min = self.cv3min
            cv3max = self.cv3max
        
        self.minima = []

        if print_output:
            print(f"Calculating gradients for FES with {self.fes.flatten().shape[0]} bins... ")
        
        if self.cvs == 1:
            m_fes = np.zeros((self.fes.shape))
            minima_count = 0
            dirs = np.zeros((self.fes.shape[0], 1))
            minima_list = []
            
            periodic=self.periodic
            
            for i in range(self.fes.shape[0]):
                a_indexes = self._get_indexes(i)
        
                okoli = np.zeros((a_indexes[:].shape))
                for ii in range(okoli.shape[0]):
                    okoli[ii] = self.fes[int(a_indexes[ii])]
                
                m_i = np.unravel_index(np.argmin(okoli), okoli.shape, order='C')
                
                dirs[i,0] = a_indexes[m_i]
                if self.fes[i] == self.fes[int(a_indexes[m_i])] and not np.all(okoli == np.min(okoli)):
                    minima_count += 1
                    m_fes[i] = minima_count
                    min_cv1 = (((i)/self.res)*(cv1max-cv1min))+cv1min
                    minima_list.append([self.fes[i],i, min_cv1])

            if print_output:
                print("Searching for the nearest local minima... ")
            iteration = 0
            while 0.0 in m_fes[:]:
                iteration += 1
                if iteration > max_iteration:
                    print("Warning: Maximum number of iterations reached when searching. ")
                    break
                for i in range(self.fes.shape[0]):
                    m_fes[i] = m_fes[int(dirs[i])]
            if iteration <= max_iteration and print_output:
                print("Done.")
            if 0.0 in m_fes[:]:
                if print_output:
                    print("Warning: some of the FES bins were not associated to any local minimum. ")

            self.minima = np.array(minima_list)
            self.m_fes = m_fes
            
            if energy_unit == "kJ/mol":
                prob = np.exp(-1000*self.fes/8.314/temp)
            elif energy_unit == "kcal_mol":
                prob = np.exp(-1000*4.184*self.fes/8.314/temp)

            m_sum_probabilities = np.zeros((self.minima[:,0].shape))
            for i in range(self.fes.shape[0]):
                if m_fes[i] > 0.0:
                    m_sum_probabilities[int(m_fes[i])-1] += prob[i]

            if energy_unit == "kJ/mol":
                m_energies = -8.314*temp*np.log(m_sum_probabilities)/1000
            elif energy_unit == "kcal_mol":
                m_energies = -8.314*temp*np.log(m_sum_probabilities)/1000/4.184

            m_energies = m_energies - np.min(m_energies)
            self.minima[:,0] = m_energies
            
        elif self.cvs == 2:
            m_fes = np.zeros((self.fes.shape))
            minima_count = 0
            dirs = np.zeros((self.fes.shape[0], self.fes.shape[1], 2))
            minima_list = []
            
            periodic=self.periodic
            
            for i in range(self.fes.shape[0]):
                for j in range(self.fes.shape[1]):
                    a_indexes = self._get_indexes(i,j)
            
                    okoli = np.zeros((a_indexes[:,:,0].shape))
                    for ii in range(okoli.shape[0]):
                        for jj in range(okoli.shape[1]):
                            okoli[ii,jj] = self.fes[int(a_indexes[ii,jj,0]),int(a_indexes[ii,jj,1])]
                    
                    m_i, m_j = np.unravel_index(np.argmin(okoli), okoli.shape, order='C')

                    if not np.all(okoli == np.min(okoli)):
                        dirs[i,j, 0] = a_indexes[m_i,m_j,0]
                        dirs[i,j, 1] = a_indexes[m_i,m_j,1]
                        if self.fes[i,j] == self.fes[int(a_indexes[m_i,m_j,0]), int(a_indexes[m_i, m_j, 1])]:
                            minima_count += 1
                            m_fes[i,j] = minima_count
                            min_cv1 = (((i)/self.res)*(cv1max-cv1min))+cv1min
                            min_cv2 = (((j)/self.res)*(cv2max-cv2min))+cv2min
                            minima_list.append([self.fes[i,j],i,j, min_cv1, min_cv2])
                    else:
                        dirs[i,j, 0] = i
                        dirs[i,j, 1] = j

            if print_output:
                print("Searching for the nearest local minima... ")
            iteration = 0
            while 0.0 in m_fes[:,:]:
                iteration += 1
                if iteration > max_iteration:
                    if print_output:
                        print("Warning: Maximum number of iterations reached when searching. ")
                    break
                for i in range(self.fes.shape[0]):
                    for j in range(self.fes.shape[1]):
                        m_fes[i,j] = m_fes[int(dirs[i,j,0]), int(dirs[i,j,1])]
            if iteration <= max_iteration:
                if print_output:
                    print("Done.")  
            if 0.0 in m_fes[:,:]:
                if print_output:
                    print("Warning: some of the FES bins were not associated to any local minimum. ")
            
            self.minima = np.array(minima_list)
            self.m_fes = m_fes
            
            if energy_unit == "kJ/mol":
                prob = np.exp(-1000*self.fes/8.314/temp)
            elif energy_unit == "kcal_mol":
                prob = np.exp(-1000*4.184*self.fes/8.314/temp)

            m_sum_probabilities = np.zeros((self.minima[:,0].shape))
            for i in range(self.fes.shape[0]):
                for j in range(self.fes.shape[1]):
                    if m_fes[i,j] > 0.0:
                        m_sum_probabilities[int(m_fes[i,j])-1] += prob[i,j]

            if energy_unit == "kJ/mol":
                m_energies = -8.314*temp*np.log(m_sum_probabilities)/1000
            elif energy_unit == "kcal_mol":
                m_energies = -8.314*temp*np.log(m_sum_probabilities)/1000/4.184

            m_energies = m_energies - np.min(m_energies)
            self.minima[:,0] = m_energies
            
        elif self.cvs == 3:
            m_fes = np.zeros((self.fes.shape))
            minima_count = 0
            dirs = np.zeros((self.fes.shape[0], self.fes.shape[1], self.fes.shape[2],3))
            minima_list = []
            
            periodic=self.periodic
            
            for i in range(self.fes.shape[0]):
                for j in range(self.fes.shape[1]):
                    for k in range(self.fes.shape[2]):
                        a_indexes = self._get_indexes(i,j,k)
                        
                        okoli = np.zeros((a_indexes[:,:,:,0].shape))
                        for ii in range(okoli.shape[0]):
                            for jj in range(okoli.shape[1]):
                                for kk in range(okoli.shape[2]):
                                    okoli[ii,jj, kk] = self.fes[int(a_indexes[ii,jj,kk,0]),int(a_indexes[ii,jj,kk,1]),int(a_indexes[ii,jj,kk,2])]
                        
                        m_i, m_j, m_k = np.unravel_index(np.argmin(okoli), okoli.shape, order='C')
                        
                        dirs[i,j,k,0] = a_indexes[m_i,m_j,m_k,0]
                        dirs[i,j,k,1] = a_indexes[m_i,m_j,m_k,1]
                        dirs[i,j,k,2] = a_indexes[m_i,m_j,m_k,2]
                        if self.fes[i,j,k] == self.fes[int(a_indexes[m_i,m_j,m_k,0]), int(a_indexes[m_i, m_j,m_k,1]), int(a_indexes[m_i, m_j,m_k,2])] and not np.all(okoli == np.min(okoli)):
                            minima_count += 1
                            m_fes[i,j,k] = minima_count
                            min_cv1 = (((i)/self.res)*(cv1max-cv1min))+cv1min
                            min_cv2 = (((j)/self.res)*(cv2max-cv2min))+cv2min
                            min_cv3 = (((k)/self.res)*(cv3max-cv3min))+cv3min
                            minima_list.append([self.fes[i,j,k],i,j,k,min_cv1,min_cv2,min_cv3])

            if print_output:
                print("Searching for the nearest local minima... ")
            iteration = 0
            while 0.0 in m_fes[:,:,:]:
                iteration += 1
                if iteration > max_iteration:
                    print("Warning: Maximum number of iterations reached when searching. ")
                    break
                for i in range(self.fes.shape[0]):
                    for j in range(self.fes.shape[1]):
                        for k in range(self.fes.shape[2]):
                            m_fes[i,j,k] = m_fes[int(dirs[i,j,k,0]), int(dirs[i,j,k,1]),int(dirs[i,j,k,2])]
            if iteration <= max_iteration:
                if print_output:
                    print("Done.")
            if 0.0 in m_fes[:,:,:]:
                if print_output:
                    print("Warning: some of the FES bins were not associated to any local minimum. ")

            self.minima = np.array(minima_list)
            self.m_fes = m_fes
            
            if energy_unit == "kJ/mol":
                prob = np.exp(-1000*self.fes/8.314/temp)
            elif energy_unit == "kcal_mol":
                prob = np.exp(-1000*4.184*self.fes/8.314/temp)

            m_sum_probabilities = np.zeros((self.minima[:,0].shape))
            for i in range(self.fes.shape[0]):
                for j in range(self.fes.shape[1]):
                    for k in range(self.fes.shape[2]):
                        if m_fes[i,j,k] > 0.0:
                            m_sum_probabilities[int(m_fes[i,j,k])-1] += prob[i,j,k]

            if energy_unit == "kJ/mol":
                m_energies = -8.314*temp*np.log(m_sum_probabilities)/1000
            elif energy_unit == "kcal_mol":
                m_energies = -8.314*temp*np.log(m_sum_probabilities)/1000/4.184

            m_energies = m_energies - np.min(m_energies)
            self.minima[:,0] = m_energies
        else:
            print("Fes object has unsupported number of CVs.")

        
        if len(self.minima.shape)>1:
            self.minima = self.minima[self.minima[:, 0].argsort()]

            # renumber minima indexes in m_fes according to their energies (and letters)
            new_m_fes = np.zeros((self.fes.shape))
            
            for m in range(self.minima.shape[0]):
                if self.cvs == 1:
                    m_old_id = self.m_fes[int(float(self.minima[m,1]))]
                if self.cvs == 2:
                    m_old_id = self.m_fes[int(float(self.minima[m,1])), int(float(self.minima[m,2]))]
                if self.cvs == 3:
                    m_old_id = self.m_fes[int(float(self.minima[m,1])), int(float(self.minima[m,2])), int(float(self.minima[m,3]))]
                new_m_fes[self.m_fes == m_old_id] = m + 1
            self.m_fes = new_m_fes

        
        letters = list(map(chr, range(65, 91)))
        for letter1 in range(65, 91):
            for letter2 in range(65, 91):
                letters.append(f"{chr(letter1)}{chr(letter2)}")
        if len(self.minima.shape) > 1:
            if self.minima.shape[0] < len(letters):
                self.minima = np.column_stack((letters[0:self.minima.shape[0]],self.minima))
            else:
                print("Error: Too many minima to assign letters. Try using keyword precise=False. ")
                return None
        elif len(self.minima.shape) == 1:
            self.minima = np.append("A", self.minima)
        
        if self.cvs == 1:
            if len(self.minima.shape)>1:
                self.minima = pd.DataFrame(np.array(self.minima), columns = ["Minimum", "free energy", "CV1bin", "CV1 - "+self.cv1_name])
            elif len(self.minima.shape) == 1:
                self.minima = pd.DataFrame([self.minima], columns = ["Minimum", "free energy", "CV1bin", "CV1 - "+self.cv1_name])
                
        elif self.cvs == 2:
            if len(self.minima.shape)>1:
                self.minima = pd.DataFrame(np.array(self.minima), columns = ["Minimum", "free energy", "CV1bin", "CV2bin", 
                                                               "CV1 - "+self.cv1_name, "CV2 - "+self.cv2_name])
            elif len(self.minima.shape) == 1:
                self.minima = pd.DataFrame([self.minima], columns = ["Minimum", "free energy", "CV1bin", "CV2bin", 
                                                               "CV1 - "+self.cv1_name, "CV2 - "+self.cv2_name])
        elif self.cvs == 3:
            if len(self.minima.shape)>1:
                self.minima = pd.DataFrame(np.array(self.minima), columns = ["Minimum", "free energy", "CV1bin", "CV2bin", "CV3bin", 
                                                               "CV1 - "+self.cv1_name, "CV2 - "+self.cv2_name,  "CV3 - "+self.cv3_name])
            elif len(self.minima.shape) == 1:
                self.minima = pd.DataFrame([self.minima], columns = ["Minimum", "free energy", "CV1bin", "CV2bin", "CV3bin", 
                                                               "CV1 - "+self.cv1_name, "CV2 - "+self.cv2_name,  "CV3 - "+self.cv3_name])


    def findminima(self, nbins=8, print_output=True):
        if int(nbins) != nbins:
            nbins = int(nbins)
            print(f"Number of bins must be an integer, it will be set to {nbins}.")
        if self.res%nbins != 0:
            print("Error: Resolution of FES must be divisible by number of bins.")
            return None
        if nbins > self.res/2:
            print("Error: Number of bins is too high.")
            return None
        bin_size = int(self.res/nbins)

        if self.cvs >= 1:
            cv1min = self.cv1min
            cv1max = self.cv1max
        if self.cvs >=2:
            cv2min = self.cv2min
            cv2max = self.cv2max 
        if self.cvs == 3:
            cv3min = self.cv3min
            cv3max = self.cv3max
        
        self.minima = []
        if self.cvs == 1:
            for bin1 in range(0,nbins):
                fes_slice = self.fes[bin1*bin_size:(bin1+1)*bin_size]
                bin_min = np.min(fes_slice)
                argmin = np.argmin(fes_slice)
                # indexes of global minimum of a bin
                bin_min_arg_cv1 = int(argmin%bin_size)
                # indexes of that minima in the original fes (indexes +1)
                min_cv1_b = int(bin_min_arg_cv1+bin1*bin_size)
                if (bin_min_arg_cv1 > 0 and bin_min_arg_cv1<(bin_size-1)):
                    min_cv1 = (((min_cv1_b)/self.res)*(cv1max-cv1min))+cv1min#min_cv1 = (((min_cv1_b+0.5)/self.res)*(cv1max-cv1min))+cv1min
                    if len(self.minima) == 0:
                        self.minima = np.array([round(bin_min, 6), int(min_cv1_b), round(min_cv1, 6)])
                    else:
                        self.minima = np.vstack((self.minima, np.array([round(bin_min, 6), int(min_cv1_b), round(min_cv1, 6)])))
                else:
                    around = []
                    min_cv1_b_low = min_cv1_b - 1
                    if min_cv1_b_low == -1:
                        if self.periodic[0]:
                            min_cv1_b_low = self.res - 1
                        else:
                            min_cv1_b_low = float("nan")

                    min_cv1_b_high = min_cv1_b + 1
                    if min_cv1_b_high == self.res:
                        if self.periodic[0]:
                            min_cv1_b_high = 0
                        else:
                            min_cv1_b_high = float("nan")

                    #1_b_low
                    if not(np.isnan(min_cv1_b_low)):
                        around.append(self.fes[min_cv1_b_low])
                    #1_b_high
                    if not(np.isnan(min_cv1_b_high)):
                        around.append(self.fes[min_cv1_b_high])
                    
                    if bin_min < np.min(around):
                        min_cv1 = (((min_cv1_b)/self.res)*(cv1max-cv1min))+cv1min#min_cv1 = (((min_cv1_b+0.5)/self.res)*(cv1max-cv1min))+cv1min
                        if len(self.minima) == 0:
                            self.minima=np.array([round(bin_min, 6), int(min_cv1_b), round(min_cv1, 6)])
                        else:
                            self.minima=np.vstack((self.minima, np.array([round(bin_min, 6), int(min_cv1_b), round(min_cv1, 6)])))
            
        elif self.cvs == 2:
            for bin1 in range(0,nbins):
                for bin2 in range(0,nbins):
                    fes_slice = self.fes[bin1*bin_size:(bin1+1)*bin_size,
                                         bin2*bin_size:(bin2+1)*bin_size]
                    bin_min = np.min(fes_slice)
                    argmin = np.argmin(fes_slice)
                    # indexes of global minimum of a bin
                    bin_min_arg = np.unravel_index(np.argmin(fes_slice), fes_slice.shape)
                    # indexes of that minima in the original fes (indexes +1)
                    min_cv1_b = int(bin_min_arg[0]+bin1*bin_size)
                    min_cv2_b = int(bin_min_arg[1]+bin2*bin_size)
                    if (bin_min_arg[0] > 0 and bin_min_arg[0]<(bin_size-1)) \
                                    and (bin_min_arg[1] > 0 and bin_min_arg[1]<(bin_size-1)):
                        min_cv1 = (((min_cv1_b)/self.res)*(cv1max-cv1min))+cv1min#min_cv1 = (((min_cv1_b+0.5)/self.res)*(cv1max-cv1min))+cv1min
                        min_cv2 = (((min_cv2_b)/self.res)*(cv2max-cv2min))+cv2min#min_cv1 = (((min_cv1_b+0.5)/self.res)*(cv1max-cv1min))+cv1min
                        if len(self.minima) == 0:
                            self.minima=np.array([round(bin_min, 6), int(min_cv1_b),\
                                                  int(min_cv2_b), round(min_cv1, 6), round(min_cv2, 6)])
                        else:
                            self.minima=np.vstack((self.minima, np.array([round(bin_min, 6), int(min_cv1_b), \
                                                                          int(min_cv2_b), round(min_cv1, 6), round(min_cv2, 6)])))
                    else:
                        around = []
                        min_cv1_b_low = min_cv1_b - 1
                        if min_cv1_b_low == -1:
                            if self.periodic[0]:
                                min_cv1_b_low = self.res - 1
                            else:
                                min_cv1_b_low = float("nan")

                        min_cv1_b_high = min_cv1_b + 1
                        if min_cv1_b_high == self.res:
                            if self.periodic[0]:
                                min_cv1_b_high = 0
                            else:
                                min_cv1_b_high = float("nan")

                        min_cv2_b_low = min_cv2_b - 1
                        if min_cv2_b_low == -1:
                            if self.periodic[0]:
                                min_cv2_b_low = self.res - 1
                            else:
                                min_cv2_b_low = float("nan")

                        min_cv2_b_high = min_cv2_b + 1
                        if min_cv2_b_high == self.res:
                            if self.periodic[0]:
                                min_cv2_b_high = 0
                            else:
                                min_cv2_b_high = float("nan")
                        #1_b_low
                        if not(np.isnan(min_cv1_b_low)):
                            if not(np.isnan(min_cv2_b_low)):
                                around.append(self.fes[min_cv1_b_low, min_cv2_b_low])
                            around.append(self.fes[min_cv1_b_low,min_cv2_b])
                            if not(np.isnan(min_cv2_b_high)):
                                around.append(self.fes[min_cv1_b_low, min_cv2_b_high])
                        #1_b
                        if not(np.isnan(min_cv2_b_low)):
                            around.append(self.fes[min_cv1_b, min_cv2_b_low])
                        if not(np.isnan(min_cv2_b_high)):
                            around.append(self.fes[min_cv1_b, min_cv2_b_high])
                        #1_b_high
                        if not(np.isnan(min_cv1_b_high)):
                            if not(np.isnan(min_cv2_b_low)):
                                around.append(self.fes[min_cv1_b_high, min_cv2_b_low])
                            around.append(self.fes[min_cv1_b_high, min_cv2_b])
                            if not(np.isnan(min_cv2_b_high)):
                                around.append(self.fes[min_cv1_b_high, min_cv2_b_high])
                        if bin_min < np.min(around):
                            min_cv1 = (((min_cv1_b)/self.res)*(cv1max-cv1min))+cv1min#min_cv1 = (((min_cv1_b+0.5)/self.res)*(cv1max-cv1min))+cv1min
                            min_cv2 = (((min_cv2_b)/self.res)*(cv2max-cv2min))+cv2min#min_cv1 = (((min_cv1_b+0.5)/self.res)*(cv1max-cv1min))+cv1min
                            if len(self.minima) == 0:
                                self.minima=np.array([round(bin_min, 6), int(min_cv1_b), int(min_cv2_b), \
                                                      round(min_cv1, 6), round(min_cv2, 6)])
                            else:
                                self.minima=np.vstack((self.minima, np.array([round(bin_min, 6), int(min_cv1_b), \
                                                                              int(min_cv2_b), round(min_cv1, 6), round(min_cv2, 6)])))
        elif self.cvs == 3:
            for bin1 in range(0,nbins):
                for bin2 in range(0,nbins):
                    for bin3 in range(0, nbins):
                        fes_slice = self.fes[bin1*bin_size:(bin1+1)*bin_size,
                                             bin2*bin_size:(bin2+1)*bin_size, 
                                             bin3*bin_size:(bin3+1)*bin_size]
                        bin_min = np.min(fes_slice)
                        argmin = np.argmin(fes_slice)
                        # indexes of global minimum of a bin
                        bin_min_arg = np.unravel_index(np.argmin(fes_slice), fes_slice.shape)
                        # indexes of that minima in the original fes (indexes +1)
                        min_cv1_b = int(bin_min_arg[0]+bin1*bin_size)
                        min_cv2_b = int(bin_min_arg[1]+bin2*bin_size)
                        min_cv3_b = int(bin_min_arg[2]+bin3*bin_size)
                        if (bin_min_arg[0] > 0 and bin_min_arg[0]<(bin_size-1)) \
                                        and (bin_min_arg[1] > 0 and bin_min_arg[1]<(bin_size-1))\
                                        and (bin_min_arg[2] > 0 and bin_min_arg[2]<(bin_size-1)):
                            min_cv1 = (((min_cv1_b)/self.res)*(cv1max-cv1min))+cv1min #min_cv1 = (((min_cv1_b+0.5)/self.res)*(cv1max-cv1min))+cv1min
                            min_cv2 = (((min_cv2_b)/self.res)*(cv2max-cv2min))+cv2min #min_cv2 = (((min_cv2_b+0.5)/self.res)*(cv2max-cv2min))+cv2min
                            min_cv3 = (((min_cv3_b)/self.res)*(cv3max-cv3min))+cv3min #min_cv3 = (((min_cv3_b+0.5)/self.res)*(cv3max-cv3min))+cv3min
                            if len(self.minima) == 0:
                                self.minima=np.array([round(bin_min, 6), int(min_cv1_b),\
                                                      int(min_cv2_b), int(min_cv3_b), round(min_cv1, 6), \
                                                      round(min_cv2, 6), round(min_cv3, 6)])
                            else:
                                self.minima=np.vstack((self.minima, np.array([round(bin_min, 6), int(min_cv1_b),\
                                                      int(min_cv2_b), int(min_cv3_b), round(min_cv1, 6), \
                                                      round(min_cv2, 6), round(min_cv3, 6)])))
                        else:
                            around = []
                            min_cv1_b_low = min_cv1_b - 1
                            if min_cv1_b_low == -1:
                                if self.periodic[0]:
                                    min_cv1_b_low = self.res - 1
                                else:
                                    min_cv1_b_low = float("nan")

                            min_cv1_b_high = min_cv1_b + 1
                            if min_cv1_b_high == self.res:
                                if self.periodic[0]:
                                    min_cv1_b_high = 0
                                else:
                                    min_cv1_b_high = float("nan")

                            min_cv2_b_low = min_cv2_b - 1
                            if min_cv2_b_low == -1:
                                if self.periodic[0]:
                                    min_cv2_b_low = self.res - 1
                                else:
                                    min_cv2_b_low = float("nan")

                            min_cv2_b_high = min_cv2_b + 1
                            if min_cv2_b_high == self.res:
                                if self.periodic[0]:
                                    min_cv2_b_high = 0
                                else:
                                    min_cv2_b_high = float("nan")
                                                       
                            min_cv3_b_low = min_cv3_b - 1
                            if min_cv3_b_low == -1:
                                if self.periodic[2]:
                                    min_cv3_b_low = self.res - 1
                                else:
                                    min_cv3_b_low = float("nan")

                            min_cv3_b_high = min_cv3_b + 1
                            if min_cv3_b_high == self.res:
                                if self.periodic[2]:
                                    min_cv3_b_high = 0
                                else:
                                    min_cv3_b_high = float("nan")

#cv3_b
                            #1_b_low
                            if not(np.isnan(min_cv1_b_low)):
                                if not(np.isnan(min_cv2_b_low)):
                                    around.append(self.fes[min_cv1_b_low,min_cv2_b_low,min_cv3_b])
                                around.append(self.fes[min_cv1_b_low,min_cv2_b,min_cv3_b])
                                if not(np.isnan(min_cv2_b_high)):
                                    around.append(self.fes[min_cv1_b_low,min_cv2_b_high,min_cv3_b])
                            #1_b
                            if not(np.isnan(min_cv2_b_low)):
                                around.append(self.fes[min_cv1_b,min_cv2_b_low,min_cv3_b])
                            if not(np.isnan(min_cv2_b_high)):
                                around.append(self.fes[min_cv1_b,min_cv2_b_high,min_cv3_b])
                            #1_b_high
                            if not(np.isnan(min_cv1_b_high)):
                                if not(np.isnan(min_cv2_b_low)):
                                    around.append(self.fes[min_cv1_b_high,min_cv2_b_low,min_cv3_b])
                                around.append(self.fes[min_cv1_b_high,min_cv2_b,min_cv3_b])
                                if not(np.isnan(min_cv2_b_high)):
                                    around.append(self.fes[min_cv1_b_high,min_cv2_b_high,min_cv3_b])
                           
                            if not(np.isnan(min_cv3_b_low)):
                            #1_b_low
                                if not(np.isnan(min_cv1_b_low)):
                                    if not(np.isnan(min_cv2_b_low)):
                                        around.append(self.fes[min_cv1_b_low,min_cv2_b_low,min_cv3_b_low])
                                    around.append(self.fes[min_cv1_b_low,min_cv2_b,min_cv3_b_low])
                                    if not(np.isnan(min_cv2_b_high)):
                                        around.append(self.fes[min_cv1_b_low,min_cv2_b_high,min_cv3_b_low])
                                #1_b
                                if not(np.isnan(min_cv2_b_low)):
                                    around.append(self.fes[min_cv1_b,min_cv2_b_low,min_cv3_b_low])
                                if not(np.isnan(min_cv2_b_high)):
                                    around.append(self.fes[min_cv1_b,min_cv2_b_high,min_cv3_b_low])
                                #1_b_high
                                if not(np.isnan(min_cv1_b_high)):
                                    if not(np.isnan(min_cv2_b_low)):
                                        around.append(self.fes[min_cv1_b_high,min_cv2_b_low,min_cv3_b_low])
                                    around.append(self.fes[min_cv1_b_high,min_cv2_b,min_cv3_b_low])
                                    if not(np.isnan(min_cv2_b_high)):
                                        around.append(self.fes[min_cv1_b_high,min_cv2_b_high,min_cv3_b_low])
                            
                            if not(np.isnan(min_cv2_b_high)):
                                #1_b_low
                                if not(np.isnan(min_cv1_b_low)):
                                    if not(np.isnan(min_cv2_b_low)):
                                        around.append(self.fes[min_cv1_b_low,min_cv2_b_low,min_cv3_b_high])
                                    around.append(self.fes[min_cv1_b_low,min_cv2_b,min_cv3_b_high])
                                    if not(np.isnan(min_cv2_b_high)):
                                        around.append(self.fes[min_cv1_b_low,min_cv2_b_high,min_cv3_b_high])
                                #1_b
                                if not(np.isnan(min_cv2_b_low)):
                                    around.append(self.fes[min_cv1_b,min_cv2_b_low,min_cv3_b_high])
                                if not(np.isnan(min_cv2_b_high)):
                                    around.append(self.fes[min_cv1_b,min_cv2_b_high,min_cv3_b_high])
                                #1_b_high
                                if not(np.isnan(min_cv1_b_high)):
                                    if not(np.isnan(min_cv2_b_low)):
                                        around.append(self.fes[min_cv1_b_high,min_cv2_b_low,min_cv3_b_high])
                                    around.append(self.fes[min_cv1_b_high,min_cv2_b,min_cv3_b_high])
                                    if not(np.isnan(min_cv2_b_high)):
                                        around.append(self.fes[min_cv1_b_high,min_cv2_b_high,min_cv3_b_high])
                            
                            if bin_min < np.min(around):
                                min_cv1 = (((min_cv1_b)/self.res)*(cv1max-cv1min))+cv1min #min_cv1 = (((min_cv1_b+0.5)/self.res)*(cv1max-cv1min))+cv1min
                                min_cv2 = (((min_cv2_b)/self.res)*(cv2max-cv2min))+cv2min #min_cv2 = (((min_cv2_b+0.5)/self.res)*(cv2max-cv2min))+cv2min
                                min_cv3 = (((min_cv3_b)/self.res)*(cv3max-cv3min))+cv3min #min_cv3 = (((min_cv3_b+0.5)/self.res)*(cv3max-cv3min))+cv3min
                                if len(self.minima) == 0:
                                    self.minima=np.array([round(bin_min, 6), int(min_cv1_b),\
                                                      int(min_cv2_b), int(min_cv3_b), round(min_cv1, 6), \
                                                      round(min_cv2, 6), round(min_cv3, 6)])
                                else:
                                    self.minima=np.vstack((self.minima, np.array([round(bin_min, 6), int(min_cv1_b),\
                                                      int(min_cv2_b), int(min_cv3_b), round(min_cv1, 6), \
                                                      round(min_cv2, 6), round(min_cv3, 6)])))
        
        else:
            print("Fes object has unsupported number of CVs.")
        
        if len(self.minima.shape)>1:
            self.minima = self.minima[self.minima[:, 0].argsort()]

        letters = list(map(chr, range(65, 91)))
        for letter1 in range(65, 91):
            for letter2 in range(65, 91):
                letters.append(f"{chr(letter1)}{chr(letter2)}")
        if len(self.minima.shape)>1:
            if self.minima.shape[1] < len(letters):
                self.minima = np.column_stack((letters[0:self.minima.shape[0]],self.minima))
            else:
                print("Error: Too many minima to assign letters.")
        elif len(self.minima.shape) == 1:
            self.minima = np.append("A", self.minima)
        
        if self.cvs == 1:
            if len(self.minima.shape)>1:
                self.minima = pd.DataFrame(np.array(self.minima), columns = ["Minimum", "free energy", "CV1bin", "CV1 - "+self.cv1_name])
            elif len(self.minima.shape) == 1:
                self.minima = pd.DataFrame([self.minima], columns = ["Minimum", "free energy", "CV1bin", "CV1 - "+self.cv1_name])
                
        elif self.cvs == 2:
            if len(self.minima.shape)>1:
                self.minima = pd.DataFrame(np.array(self.minima), columns = ["Minimum", "free energy", "CV1bin", "CV2bin", 
                                                               "CV1 - "+self.cv1_name, "CV2 - "+self.cv2_name])
            elif len(self.minima.shape) == 1:
                self.minima = pd.DataFrame([self.minima], columns = ["Minimum", "free energy", "CV1bin", "CV2bin", 
                                                               "CV1 - "+self.cv1_name, "CV2 - "+self.cv2_name])
        elif self.cvs == 3:
            if len(self.minima.shape)>1:
                self.minima = pd.DataFrame(np.array(self.minima), columns = ["Minimum", "free energy", "CV1bin", "CV2bin", "CV3bin", 
                                                               "CV1 - "+self.cv1_name, "CV2 - "+self.cv2_name,  "CV3 - "+self.cv3_name])
            elif len(self.minima.shape) == 1:
                self.minima = pd.DataFrame([self.minima], columns = ["Minimum", "free energy", "CV1bin", "CV2bin", "CV3bin", 
                                                               "CV1 - "+self.cv1_name, "CV2 - "+self.cv2_name,  "CV3 - "+self.cv3_name])
        

    def plot(self, png_name=None, contours=True, contours_spacing=0.0, aspect = 1.0, cmap = "RdYlBu_r", 
                 energy_unit="kJ/mol", xlabel=None, ylabel=None, zlabel=None, label_size=12, clabel_size = 12, image_size=None, image_size_unit="in", dpi=100, color=None, vmin = 0, vmax = None, opacity=0.2, levels=None, show_points=True, point_size=4.0, title = None, off_screen = False, xlim=[None, None], ylim=[None, None], return_fig=False):
        """
        The same function as for visualizing Fes objects, but this time 
        with the positions of local minima shown as letters on the graph.
        
        ```python
        minima.plot(png_name="minima.png")
        ```
        
        Parameters:
        
        * png_name = String. If this parameter is supplied, the picture of FES will be saved under this name to the current working directory.
        
        * contours (default=True) = whether contours should be shown on 2D FES
        
        * contours_spacing (default=0.0) = when a positive number is set, it will be used as spacing for contours on 2D FES. 
                Otherwise, if contours=True, there will be five equally spaced contour levels.
        
        * aspect (default = 1.0) = aspect ratio of the graph. Works with 1D and 2D FES. 
        
        * cmap (default = "RdYlBu_r") = Matplotlib colormap used to color 2D or 3D FES
        
        * energy_unit (default="kJ/mol") = String, used in description of colorbar

        * xlim, ylim (default = [None, None] for both) = list of two values specifying the range of x and y axes respectively. 
        None means that Matplotlib will choose appropriate range, so this keyword is useful when you need to overwrite it. 
        
        * xlabel, ylabel, zlabel = Strings, if provided, they will be used as labels for the graphs
        
        * label_size, clabel_size (default = 12) = size of text in labels or contours labels, respectively
        
        * image_size (default = [9,6]) = List of the width and height of the picture

        * image_size_unit (default = "in") = Units for width and height of the picture, accepts "in" as inches, "cm", "mm" or "px" as pixels.

        * dpi (default = 100) = DPI of the resulting image. 
        
        * color = string = name of color in matplotlib, if set, the color will be used for the letters. 
                If not set, the color should be automatically either black or white, 
                depending on what will be better visible on given place on FES with given colormap (for 2D FES).
        
        * vmin (default=0) = real number, lower bound for the colormap on 2D FES
        
        * vmax = real number, upper bound for the colormap on 2D FES
        
        * opacity (default=0.2) = number between 0 and 1, is the opacity of isosurfaces of 3D FES
        
        * levels = Here you can specify list of free energy values for isosurfaces on 3D FES. 
                If not provided, default values from contours parameters will be used instead. 
        
        * show_points (default=True) = boolean, tells if points should be visualized too, instead of just the letters. Only on 3D FES. 
        
        * point_size (default=4.0) = float, sets the size of points if show_points=True

        * title = optional, string that defines the title of the graph

        * offscreen (default = False) = for 3D FES only, grapf will not be shown after creation - used internally when making animations

        * return_fig (default=False) = whether the method should return the Matplotlib.Pyplot.figure object for further use. 
        In the case of plotting 3D FES, it returns Pyvista.plotter object instead. 
        """
        
        if vmax == None:
            vmax = np.max(self.fes)+0.01 # if the addition is smaller than 0.01, the 3d plot stops working. 
            
        if contours_spacing == 0.0:
            contours_spacing = (vmax-vmin)/5.0
        
        cmap = cm.get_cmap(cmap)
        
        cmap.set_over("white")
        cmap.set_under("white")
        
        color_set = True
        if color == None:
            color_set = False
        
        if self.cvs >= 1:
            cv1min = self.cv1min
            cv1max = self.cv1max 
        if self.cvs >=2:
            cv2min = self.cv2min
            cv2max = self.cv2max
        if self.cvs == 3:
            cv3min = self.cv3min
            cv3max = self.cv3max

        if image_size == None:
            image_size = [9,6]
        
        if image_size_unit == "cm":
            image_size[0] /= 2.54
            image_size[1] /= 2.54
        elif image_size_unit == "mm":
            image_size[0] /= 25.4
            image_size[1] /= 25.4
        elif image_size_unit == "px":
            image_size[0] /= dpi
            image_size[1] /= dpi
        elif image_size_unit != "in":
            print(f"Warning: unknown image_size_unit value: {image_size_unit}. Using inches instead. ")
            
        if self.cvs == 1:
            fig = plt.figure(figsize=(image_size[0],image_size[1]), dpi=dpi)
            X = np.linspace(cv1min, cv1max, self.res)
            plt.plot(X, self.fes)
            
            if not color_set:
                color = "black"
            
            fesrange = np.max(self.fes) - np.min(self.fes)
            
            if self.minima.shape[0] == 1:
                plt.text(float(self.minima.iloc[0,3]), float(self.minima.iloc[0,1])+fesrange*0.05, self.minima.iloc[0,0],
                             fontsize=label_size, horizontalalignment='center',
                             verticalalignment='bottom', c=color)
            elif self.minima.shape[0] > 1:
                for m in range(len(self.minima.iloc[:,0])):
                    plt.text(float(self.minima.iloc[m,3]), float(self.minima.iloc[m,1])+fesrange*0.05, self.minima.iloc[m,0],
                             fontsize=label_size, horizontalalignment='center',
                             verticalalignment='bottom', c=color)
            
            if xlabel == None:
                plt.xlabel(f'CV1 - {self.cv1_name}', size=label_size)
            else:
                plt.xlabel(xlabel, size=label_size)
            if ylabel == None:
                plt.ylabel(f'free energy ({energy_unit})', size=label_size)
            else:
                plt.ylabel(ylabel, size=label_size)
            
            if title != None:
                plt.title(title, size=label_size)
    
            plt.xticks(fontsize=label_size)
            plt.yticks(fontsize=label_size)

            ax = plt.gca()
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
        
            if png_name != None:
                plt.savefig(png_name, bbox_inches = 'tight')

            if return_fig:
                return fig
                
            
        elif self.cvs == 2:
            fig = plt.figure(figsize=(image_size[0],image_size[1]), dpi=dpi)
            plt.imshow(np.rot90(self.fes, axes=(0,1)), cmap=cmap, interpolation='nearest', 
                       extent=[cv1min, cv1max, cv2min, cv2max], 
                       aspect = (((cv1max-cv1min)/(cv2max-cv2min))/(aspect)),
                       vmin = vmin, vmax = vmax)
            cbar = plt.colorbar()
            cbar.ax.tick_params(labelsize=label_size) 
            cbar.set_label(energy_unit, size=label_size)

            if contours:
                if levels != None:
                    cont = plt.contour(np.rot90(self.fes, axes=(0,1)), 
                             levels = levels, 
                             extent=[cv1min, cv1max, cv2max, cv2min], 
                             colors = "k", linestyles = "solid")
                    plt.clabel(cont, levels = levels, fontsize=clabel_size)
                else:
                    cont = plt.contour(np.rot90(self.fes, axes=(0,1)), 
                             levels = np.arange(vmin, (vmax - 0.01), contours_spacing), 
                             extent=[cv1min, cv1max, cv2max, cv2min], 
                             colors = "k", linestyles = "solid")
                    plt.clabel(cont, levels = np.arange(vmin, (vmax - 0.01), contours_spacing), fontsize=clabel_size)
            
            if self.minima.shape[0] == 1:
                background = cmap((float(self.minima.iloc[0,1])-vmin)/(vmax-vmin))
                luma = background[0]*0.2126+background[1]*0.7152+background[3]*0.0722
                if luma > 0.6 and not color_set:
                    color = "black"
                elif luma <= 0.6 and not color_set:
                    color="white"
                plt.text(float(self.minima.iloc[0,4])+0.5*(cv1max-cv1min)/self.res, float(self.minima.iloc[0,5])+0.5*(cv2max-cv2min)/self.res, self.minima.iloc[0,0], fontsize=label_size, horizontalalignment='center',verticalalignment='center', c=color)
            elif self.minima.shape[0] > 1:
                for m in range(len(self.minima.iloc[:,0])):
                    background = cmap((float(self.minima.iloc[m,1])-vmin)/(vmax-vmin))
                    luma = background[0]*0.2126+background[1]*0.7152+background[3]*0.0722
                    if luma > 0.6 and not color_set:
                        color = "black"
                    elif luma <= 0.6 and not color_set:
                        color="white"
                    plt.text(float(self.minima.iloc[m,4])+0.5*(cv1max-cv1min)/self.res, float(self.minima.iloc[m,5])+0.5*(cv2max-cv2min)/self.res, self.minima.iloc[m,0], fontsize=label_size, horizontalalignment='center', verticalalignment='center', c=color)
            
            
            if xlabel == None:
                plt.xlabel(f'CV1 - {self.cv1_name}', size=label_size)
            else:
                plt.xlabel(xlabel, size=label_size)
            if ylabel == None:
                plt.ylabel(f'CV2 - {self.cv2_name}', size=label_size)
            else:
                plt.ylabel(ylabel, size=label_size)
            
            if title != None:
                plt.title(title, size=label_size)
    
            plt.xticks(fontsize=label_size)
            plt.yticks(fontsize=label_size)

            ax = plt.gca()
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            ax.tick_params(axis='x', labelsize=label_size)
            ax.tick_params(axis='y', labelsize=label_size)
        
            if png_name != None:
                plt.savefig(png_name, bbox_inches = 'tight')

            if return_fig:
                return fig
            
            
        elif self.cvs == 3:
            if xlabel == None:
                xlabel = "CV1 - " + self.cv1_name
            if ylabel == None:
                ylabel = "CV2 - " + self.cv2_name
            if zlabel == None:
                zlabel = "CV3 - " + self.cv3_name
            
            min_ar = self.minima.iloc[:,5:8].values
            min_ar = min_ar.astype(np.float32)
            min_pv = pv.PolyData(min_ar)
            
            
            grid = pv.ImageData(
                dimensions=(self.res, self.res, self.res),
                spacing=((cv1max-cv1min)/self.res,(cv2max-cv2min)/self.res,(cv3max-cv3min)/self.res),
                origin=(cv1min, cv2min, cv3min)
            )
            grid["vol"] = self.fes.ravel(order="F")
            if levels == None:
                contours = grid.contour(np.arange(0, (vmax - 0.1), contours_spacing))
            else:
                contours = grid.contour(levels)
            fescolors = []
            for i in range(contours.points.shape[0]):
                fescolors.append(self.fes[int((contours.points[i,0]-cv1min)*self.res/(cv1max-cv1min)),
                                          int((contours.points[i,1]-cv2min)*self.res/(cv2max-cv2min)),
                                          int((contours.points[i,2]-cv3min)*self.res/(cv3max-cv3min))])
            #%% Visualization
            bounds = [cv1min, cv1max, cv2min, cv2max, cv3min, cv3max]
            pv.set_plot_theme('document')
            p = pv.Plotter()
            p.add_mesh(contours, scalars=fescolors, opacity=opacity, cmap=cmap, show_scalar_bar=False, interpolate_before_map=True)
            p.show_bounds(bounds=bounds, xtitle=xlabel, ytitle=ylabel, ztitle=zlabel, grid=True, location="outer")
            p.add_point_labels(min_pv, self.minima.iloc[:,0], 
                   show_points=show_points, always_visible = True, 
                   point_color="black", point_size=point_size, 
                   font_size=label_size, shape=None)
            
            if title != None:
                plt.title(title)
            if not off_screen:
                p.show()

            if png_name != None:
                p.save_graphic(png_name)

            if return_fig:
                return p
            

    def make_gif(self, gif_name=None, cmap = "RdYlBu_r", energy_unit="kJ/mol", 
                 xlabel=None, ylabel=None, zlabel=None, label_size=12, image_size=[10,7], 
                  opacity=0.2, levels=None, show_points=True, point_size=4.0, frames=64):
        """
        Function that generates animation of 3D FES showing different isosurfaces.
        
        ```python
        fes.make_gif(gif_name="FES.gif")
        ```
        
        Parameters:
        
        * gif_name (default="minima.gif") = String. Name of the gif that will be saved in the working directory.
        
        * cmap (default = "RdYlBu_r") = Matplotlib colormap used to color the 3D FES
        
        * xlabel, ylabel, zlabel = Strings, if provided, they will be used as labels for the graph
        
        * labelsize (default = 12) = size of text in labels
        
        * image_size (default = [10,7]) = List of the width and height of the picture
        
        * opacity (default = 0.2) = number between 0 and 1, is the opacity of isosurfaces of 3D FES
        
        * levels = Here you can specify list of free energy values for isosurfaces on 3D FES. 
                If not provided, default values from contours parameters will be used instead. 
        
        * frames (default = 64) = Number of frames the animation will be made of. 
        """
        if self.cvs == 3:
            values = np.linspace(np.min(self.fes)+0.1, np.max(self.fes), num=frames)
            grid = pv.ImageData(
                dimensions=(self.res, self.res, self.res),
                spacing=((self.cv1max-self.cv1min)/self.res,(self.cv2max-self.cv2min)/self.res,(self.cv3max-self.cv3min)/self.res),
                origin=(self.cv1min, self.cv2min, self.cv3min),
            )
            grid["vol"] = self.fes.ravel(order="F")
            surface = grid.contour(values[:1])
            surfaces = [grid.contour([v]) for v in values]
            surface = surfaces[0].copy()
            
            pv.set_plot_theme('document')
            plotter = pv.Plotter(off_screen=True)
            # Open a movie file
            if gif_name == None:
                gif_name = "minima_animation.gif"
            plotter.open_gif(gif_name)

            # Add initial mesh
            plotter.add_mesh(
                surface,
                opacity=0.3,
                clim=grid.get_data_range(),
                show_scalar_bar=False,
                cmap="RdYlBu_r"
            )
            plotter.add_mesh(grid.outline_corners(), color="k")
            
            if xlabel == None and ylabel == None and zlabel == None:
                plotter.show_grid(xtitle=f"CV1 - {self.cv1_name}", ytitle=f"CV2 - {self.cv2_name}", ztitle=f"CV3 - {self.cv3_name}")
            else:
                plotter.show_grid(xtitle=xlabel, ytitle=ylabel, ztitle=zlabel)
            if show_points:
                min_ar = self.minima.iloc[:,5:8].values
                min_ar = min_ar.astype(np.float32)
                min_pv = pv.PolyData(min_ar)
                plotter.add_point_labels(min_pv, self.minima.iloc[:,0], 
                               show_points=True, always_visible = True, 
                               pickable = True, point_color="black", 
                               point_size=4, font_size=label_size, shape=None)
            plotter.set_background('white')
            text = plotter.add_text(f"{values[0]:.2f}+kJ/mol", position='lower_right', font_size=label_size)
            plotter.show(auto_close=False)
            
            # Run through each frame
            for surf in range(len(surfaces)):
                surface.copy_from(surfaces[surf])
                plotter.remove_actor(text)
                text = plotter.add_text(f"{values[surf]:.2f} {energy_unit}", position='lower_right', font_size=label_size)
                plotter.write_frame()  # Write this frame
            # Run through backwards
            for surf in range(len(surfaces)-1,0,-1):
                surface.copy_from(surfaces[surf])
                plotter.remove_actor(text)
                text = plotter.add_text(f"{values[surf]:.2f} {energy_unit}", position='lower_right', font_size=12)
                plotter.write_frame()  # Write this frame

            # Be sure to close the plotter when finished
            plotter.close()
        else:
            print("Error: this method is only available for FES with 3 CVs.")


class FEProfile:
    """
    Free energy profile is a visualization of differences between local 
    minima points during metadynamics simulation. If the values seem 
    to converge to a mean value of the difference, it suggests, 
    but not fully prooves, that the obtained FES did converge to the correct shape.
    
    Command:
    ```python
    fep = metadynminer.FEProfile(minima, hillsfile)
    ```
    
    Parameters:
    
    * minima = metadynminer.Minima object
    
    * hillsfile = metadynminer.Hills object
    
    """
    
    __pdoc__["FEProfile.makefeprofile"] = False
    
    def __init__(self, minima, hills):
        self.cvs = minima.cvs
        self.res = minima.res
        self.minima = minima.minima
        self.periodic = minima.periodic
        self.heights = hills.get_heights()
        
        if self.cvs >= 1:
            self.cv1_name = minima.cv1_name
            self.cv1min = minima.cv1min
            self.cv1max = minima.cv1max
            self.cv1 = hills.get_cv1()
            self.s1 = hills.get_sigma1()
            self.cv1per = hills.get_cv1per()
            self.cv1_fes_range = minima.cv1_fes_range
        if self.cvs >= 2:
            self.cv2min = minima.cv2min
            self.cv2max = minima.cv2max
            self.cv2_name = minima.cv2_name
            self.cv2 = hills.get_cv2()
            self.s2 = hills.get_sigma2()
            self.cv2per = hills.get_cv2per()
            self.cv2_fes_range = minima.cv2_fes_range
        if self.cvs == 3:
            self.cv3min = minima.cv3min
            self.cv3max = minima.cv3max
            self.cv3_name = minima.cv3_name
            self.cv3 = hills.get_cv3()
            self.s3 = hills.get_sigma3()
            self.cv3per = hills.get_cv3per()
            self.cv3_fes_range = minima.cv3_fes_range
        
        if len(minima.minima.shape)>1:
            self.makefeprofile(hills)
        else: 
            print("There is only one local minimum on the free energy surface.")
    
    def makefeprofile(self, hills):
        """
        Internal method to calculate free energy profile. 
        """
        hillslenght = len(hills.get_cv1())
        
        if hillslenght < 256:
            profilelenght = hillslenght
            scantimes = np.array(range(hillslenght))
        else:
            profilelenght = 256
            scantimes = np.array(((hillslenght/(profilelenght))*np.array((range(1,profilelenght+1)))))
            scantimes[0:-1] -= 1
            #scantimes[-1] += 1
            scantimes = scantimes.astype(int)
        number_of_minima = self.minima.shape[0]
        
        self.feprofile = np.zeros((self.minima.Minimum.shape[0]+1))

        if self.cvs >= 1:
            cv1min = self.cv1min
            cv1max = self.cv1max         
        if self.cvs >= 2:
            cv2min = self.cv2min
            cv2max = self.cv2max           
        if self.cvs >= 3:     
            cv3min = self.cv3min
            cv3max = self.cv3max
        
        if self.cvs == 1:
            
            fes = np.zeros((self.res))
            
            lasttime = 0
            for time in scantimes:
                if time == scantimes[-1]:
                    time += 1
                for m in range(number_of_minima):#self.minima.iloc[:,3]
                    dist_cv1 = self.cv1[lasttime:time]-float(self.minima.iloc[m,3])
                    if self.periodic[0]:
                        dist_cv1[dist_cv1<-0.5*self.cv1_fes_range] += self.cv1_fes_range
                        dist_cv1[dist_cv1>+0.5*self.cv1_fes_range] -= self.cv1_fes_range

                    dp2 = dist_cv1**2/(2*self.s1[lasttime:time]**2)
                    tmp = np.zeros(dp2.shape)
                    heights = self.heights[lasttime:time]
                    tmp[dp2<6.25] = heights[dp2<6.25] * (np.exp(-dp2[dp2<6.25]) * 1.00193418799744762399 - 0.00193418799744762399)
                    #fes[int(float(self.minima.iloc[x,2]))] -= tmp.sum()
                    fes[int(float(self.minima.iloc[m,2]))] -= tmp.sum()
                #fes = fes - np.min(fes)
                profileline = [time-1]
                for m in range(number_of_minima):
                    profileline.append(fes[int(float(self.minima.iloc[m,2]))]-\
                                       fes[int(float(self.minima.iloc[0,2]))])
                self.feprofile = np.vstack([self.feprofile, profileline])

                lasttime = time
            
        elif self.cvs == 2:
            
            fes = np.zeros((self.res, self.res))
            
            lasttime = 0
            line = 0
            for time in scantimes:
                if time == scantimes[-1]:
                    time += 1
                for m in range(number_of_minima):
                    dist_cv1 = self.cv1[lasttime:time]-float(self.minima.iloc[m,4])
                    if self.periodic[0]:
                        dist_cv1[dist_cv1<-0.5*self.cv1_fes_range] += self.cv1_fes_range
                        dist_cv1[dist_cv1>+0.5*self.cv1_fes_range] -= self.cv1_fes_range
                    
                    dist_cv2 = self.cv2[lasttime:time]-float(self.minima.iloc[m,5])
                    if self.periodic[1]:
                        dist_cv2[dist_cv2<-0.5*self.cv2_fes_range] += self.cv2_fes_range
                        dist_cv2[dist_cv2>+0.5*self.cv2_fes_range] -= self.cv2_fes_range
                
                    dp2 = dist_cv1**2/(2*self.s1[lasttime:time]**2) + dist_cv2**2/(2*self.s2[lasttime:time]**2)
                    tmp = np.zeros(self.cv1[lasttime:time].shape)
                    heights = self.heights[lasttime:time]
                    tmp[dp2<6.25] = heights[dp2<6.25] * (np.exp(-dp2[dp2<6.25]) * 1.00193418799744762399 - 0.00193418799744762399)
                    fes[int(float(self.minima.iloc[m,2])),int(float(self.minima.iloc[m,3]))] -= tmp.sum()
                
                # save profile
                profileline = [time-1]
                for m in range(number_of_minima):
                    profileline.append(fes[int(float(self.minima.iloc[m,2])),int(float(self.minima.iloc[m,3]))]-\
                                       fes[int(float(self.minima.iloc[0,2])),int(float(self.minima.iloc[0,3]))])
                self.feprofile = np.vstack([self.feprofile, profileline])

                lasttime = time
        
        elif self.cvs == 3:
            
            fes = np.zeros((self.res, self.res, self.res))
            
            lasttime = 0
            for time in scantimes:
                if time == scantimes[-1]:
                    time += 1
                for m in range(number_of_minima):
                    dist_cv1 = self.cv1[lasttime:time]-float(self.minima.iloc[m,5])
                    if self.periodic[0]:
                        dist_cv1[dist_cv1<-0.5*self.cv1_fes_range] += self.cv1_fes_range
                        dist_cv1[dist_cv1>+0.5*self.cv1_fes_range] -= self.cv1_fes_range
                    
                    dist_cv2 = self.cv2[lasttime:time]-float(self.minima.iloc[m,6])
                    if self.periodic[1]:
                        dist_cv2[dist_cv2<-0.5*self.cv2_fes_range] += self.cv2_fes_range
                        dist_cv2[dist_cv2>+0.5*self.cv2_fes_range] -= self.cv2_fes_range
                    
                    dist_cv3 = self.cv3[lasttime:time]-float(self.minima.iloc[m,7])
                    if self.periodic[2]:
                        dist_cv3[dist_cv3<-0.5*self.cv3_fes_range] += self.cv3_fes_range
                        dist_cv3[dist_cv3>+0.5*self.cv3_fes_range] -= self.cv3_fes_range
            
                    dp2 = (dist_cv1**2/(2*self.s1[lasttime:time]**2) + 
                           dist_cv2**2/(2*self.s2[lasttime:time]**2) + 
                           dist_cv3**2/(2*self.s3[lasttime:time]**2))
                    tmp = np.zeros(dp2.shape)
                    heights = self.heights[lasttime:time]
                    tmp[dp2<6.25] = (heights[dp2<6.25] * (np.exp(-dp2[dp2<6.25]) * 1.00193418799744762399 - 0.00193418799744762399))
                    fes[int(float(self.minima.iloc[m,2])),
                        int(float(self.minima.iloc[m,3])),
                        int(float(self.minima.iloc[m,4]))] -= tmp.sum()
                
                # save profile
                profileline = [time-1]
                for m in range(number_of_minima):
                    profileline.append(fes[int(float(self.minima.iloc[m,2])),
                                           int(float(self.minima.iloc[m,3])),
                                           int(float(self.minima.iloc[m,4]))]-\
                                       fes[int(float(self.minima.iloc[0,2])),
                                           int(float(self.minima.iloc[0,3])),
                                           int(float(self.minima.iloc[0,4]))])
                self.feprofile = np.vstack([self.feprofile, profileline])

                lasttime = time
            
        else:
            print("Fes object doesn't have supported number of CVs.")
    
    def plot(self, png_name=None, image_size=None, image_size_unit="in", dpi=100, tu = "ps", xlabel=None, ylabel=None, label_size=12, cmap="RdYlBu_r", legend=True, xlim=[None, None], ylim=[None, None], title=None, return_fig=False):
        """
        Visualization function for free energy profiles. 
        
        ```python
        fep.plot(png_name="FEProfile.png")
        ```
        
        
        Parameters:
        
        
        * png_name (default=None) = name for image file to save the plot to
        
        * image_size (default = [9,6]) = list of the width and height of the picture
        
        * image_size_unit (default = "in") = Units for width and height of the picture, accepts "in" as inches, "cm", "mm" or "px" as pixels. 

        * dpi (default = 100) = DPI of the resulting image. 
        
        * tu (default = "ps") = string, time unit to be shown on x axis. Available options: "s", "ms", "us", "ns", "ps", "fs"
        
        * xlabel (default="time (ps)")
        
        * ylabel (default="free energy difference (kJ/mol)") 
        
        * label_size (default=12) = size of labels
        
        * cmap (default="RdYlBu_r") = matplotlib colormap used for coloring the line of the minima

        * legend (default=True) = whether there should be a matplotlib's legend in the graph

        * xlim, ylim (default = [None, None] for both) = list of two values specifying the range of x and y axes respectively. 
        None means that Matplotlib will choose appropriate range, so this keyword is useful when you need to overwrite it. 

        * title = optional, string that defines the title of the graph

        * return_fig (default=False) = whether the method should return the Matplotlib.Pyplot.figure object for further use. 
        
        """
        
        tu = TU(tu)

        if image_size == None:
            image_size = [9,6]
        
        if image_size_unit == "cm":
            image_size[0] /= 2.54
            image_size[1] /= 2.54
        elif image_size_unit == "mm":
            image_size[0] /= 25.4
            image_size[1] /= 25.4
        elif image_size_unit == "px":
            image_size[0] /= dpi
            image_size[1] /= dpi
        elif image_size_unit != "in":
            print(f"Warning: unknown image_size_unit value: {image_size_unit}. Using inches instead. ")
        
        fig = plt.figure(figsize=(image_size[0],image_size[1]), dpi=dpi)
        
        cmap=cm.get_cmap(cmap)
        
        #colors = cm.RdYlBu_r((self.minima.iloc[:,1].to_numpy()).astype(float)/\
        #                (np.max(self.minima.iloc[:,1].to_numpy().astype(float))))
        colors = cmap(np.linspace(0,1,self.minima.shape[0]))
        for m in range(self.minima.shape[0]):
            plt.plot(tu.intu(self.feprofile[:,0]), self.feprofile[:,m+1], color=colors[m])

        if xlabel == None:
            plt.xlabel(f'time ({tu.name})', size=label_size)
        else:
            plt.xlabel(xlabel, size=label_size)
        if ylabel == None:
            plt.ylabel('free energy difference (kJ/mol)', size=label_size)
        else:
            plt.ylabel(ylabel, size=label_size)
        if legend:
            plt.legend(self.minima.iloc[:,0], loc="lower right")
        
        if title != None:
            plt.title(title, size=label_size)
    
        plt.xticks(fontsize=label_size)
        plt.yticks(fontsize=label_size)
        
        ax = plt.gca()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        if png_name != None:
            plt.savefig(png_name, bbox_inches = 'tight')

        if return_fig:
            return fig







