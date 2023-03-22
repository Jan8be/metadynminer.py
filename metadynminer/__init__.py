name = "metadynminer"
try:
    import numpy as np
except:
    print("Error while loading numpy")
    exit(0)
try:
    from matplotlib import pyplot as plt
except:
    print("Error while loading matplotlib pyplot")
    exit(0)
try:
    from matplotlib import cm
except:
    print("Error while loading matplotlib cm")
    exit(0)
try:
    import pandas as pd
except:
    print("Error while loading pandas")
    exit(0)
try:
    from mpl_toolkits import mplot3d
except:
    print("Error while loading mpl_toolkits")
    exit(0)


class Hills:   
    def __init__(self, name="HILLS", encoding="utf8", ignoretime=True, periodic=[False, False], 
                 cv1per=[-np.pi, np.pi],cv2per=[-np.pi, np.pi],cv3per=[-np.pi, np.pi]):
        self.read(name, encoding, ignoretime, periodic, 
                 cv1per=[-np.pi, np.pi],cv2per=[-np.pi, np.pi],cv3per=[-np.pi, np.pi])
        self.hillsfilename = name
    
    def read(self, name="HILLS", encoding="utf8", ignoretime=True, periodic=None, 
                 cv1per=[-np.pi, np.pi],cv2per=[-np.pi, np.pi],cv3per=[-np.pi, np.pi]):
        with open(name, 'r', encoding=encoding) as hillsfile:
            lines = hillsfile.readlines()
        columns = lines[0].split() 
        number_of_columns_head = len(columns) - 2
        
        if number_of_columns_head == 5:
            self.cvs = 1
            self.cv1_name = lines[0].split()[3]
            if periodic == None:
                periodic = list(False)
            
            self.periodic = list(periodic[0:1])
            self.cv1per = cv1per

        elif number_of_columns_head == 7:
            self.cvs = 2
            self.cv1_name = lines[0].split()[3]
            self.cv2_name = lines[0].split()[4]
            if len(periodic) == 2:
                self.periodic = periodic[0:2]
                self.cv1per = cv1per
                self.cv2per = cv2per
            else:
                print(f"Argument 'periodic' has wrong number of parameters({len(periodic)})")
        elif number_of_columns_head == 9:
            self.cvs = 3
            self.cv1_name = lines[0].split()[3]
            self.cv2_name = lines[0].split()[4]
            self.cv2_name = lines[0].split()[5]
            
            if len(periodic) == 3:
                self.periodic = periodic[0:3]
                self.cv1per = cv1per
                self.cv2per = cv2per
                self.cv3per = cv3per
            else:
                print(f"Argument 'periodic' has wrong number of parameters({len(periodic)})")
        else:
            print("Unexpected number of columns in provided HILLS file.")
            
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
                else:
                    self.hills.append(lines[line].split())
                if ignoretime:
                    self.hills[(t-1)][0] = t
        
        self.hills = np.array(self.hills, dtype=np.double)
                
        if self.cvs == 1:
            self.cv1 = self.hills[:,1]
            self.sigma1 = self.hills[:,2]
            self.heights = self.hills[:,3]
        elif self.cvs == 2:
            self.cv1 = self.hills[:,1]
            self.cv2 = self.hills[:,2]
            self.sigma1 = self.hills[:,3]
            self.sigma2 = self.hills[:,4]
            self.heights = self.hills[:,5]
        elif self.cvs == 3:
            self.cv1 = self.hills[:,1]
            self.cv2 = self.hills[:,2]
            self.cv3 = self.hills[:,3]
            self.sigma1 = self.hills[:,4]
            self.sigma2 = self.hills[:,5]
            self.sigma3 = self.hills[:,6]
            self.heights = self.hills[:,7]
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


class Fes: 
    def __init__(self, hills=None, resolution=256, original = False, calculate_new_fes = True):
        if hills != None:
            self.hills = hills
            self.cvs = hills.get_number_of_cvs()
            self.heights = hills.get_heights()
            self.periodic = hills.get_periodic()
            
            if self.cvs >= 1:
                self.cv1 = hills.get_cv1()
                self.s1 = hills.get_sigma1()

                self.cv1min = np.min(self.cv1) - 1e-8
                self.cv1max = np.max(self.cv1) + 1e-8

                self.cv1_name = hills.get_cv1_name()
                self.cv1per = hills.get_cv1per()
                
                if not original:
                    if ((np.max(self.s1)/np.min(self.s1))>1.00000001):
                        print("""Error: Bias sum algorithm only works for hills files 
                        in which all hills have the same width. 
                        For this file, you need the slower but exact, algorithm, to do that, 
                        set the argument 'original' to True.""")

            if self.cvs >= 2:
                self.cv2 = hills.get_cv2()
                self.s2 = hills.get_sigma2()

                self.cv2min = np.min(self.cv2) - 1e-8
                self.cv2max = np.max(self.cv2) + 1e-8

                self.cv2_name = hills.get_cv2_name()
                self.cv2per = hills.get_cv2per()
                
                if not original:
                    if ((np.max(self.s2)/np.min(self.s2))>1.00000001):
                        print("""Error: Bias sum algorithm only works for hills files 
                        in which all hills have the same width. 
                        For this file, you need the slower but exact, algorithm, to do that, 
                        set the argument 'original' to True.""")

            elif self.cvs == 3:
                self.cv3 = hills.get_cv3()
                self.s3 = hills.get_sigma3()

                self.cv3min = np.min(self.cv3) - 1e-8
                self.cv3max = np.max(self.cv3) + 1e-8

                self.cv3_name = hills.get_cv3_name()
                self.cv3per = hills.get_cv3per()
                
                if not original:
                    if ((np.max(self.s3)/np.min(self.s3))>1.00000001):
                        print("""Error: Bias sum algorithm only works for hills files 
                        in which all hills have the same sigma. 
                        For this file, you need the slower but exact, algorithm, to do that, 
                        set the argument 'original' to True.""")
            
            if not original:
                if calculate_new_fes:
                    try: self.fes
                    except AttributeError: 
                        if original:
                            ...
                        else:
                            self.makefes(self.hills, resolution)
            else:
                ...#self.makefes2(self.hills, resolution)
        
    def makefes(self, hills, resolution=256):
        self.res = resolution
        #if self.res % 2 == 0:
        #    self.res += 1
                
        if self.cvs == 1:
            cv1min = self.cv1min
            cv1max = self.cv1max
            if not self.periodic[0]:
                print("increase non periodic range")
                cv1range = self.cv1max-self.cv1min
                cv1min -= cv1range*0.15          
                cv1max += cv1range*0.15
                
            cv1bin = np.ceil((self.cv1-cv1min)*self.res/(cv1max-cv1min))
            cv1bin = cv1bin.astype(int)
            s1res = (self.s1[0]*self.res)/(cv1max - cv1min)
            
            gauss_res = 8*s1res
            gauss_res = int(gauss_res)
            if gauss_res%2 == 0:
                gauss_res += 1
            
            x = np.arange((-(gauss_res-1)/2), ((gauss_res-1)/2)+1)
            gauss = -np.exp(-x*x/2.0/s1res/s1res)
            
            fes = np.zeros((self.res))
            
            for line in range(len(cv1bin)):
                if (line) % 500 == 0:
                    print(f"Constructing free energy surface: {((line+1)/len(cv1bin)):.1%} finished", end="\r")
                
                fes_center = int((self.res-1)/2)
                gauss_center_to_end = int((gauss_res-1)/2)
                #print(f"\ng_res: {gauss_res}, gauss_center_to_end: {gauss_center_to_end}")
                
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
            print("\n")
            fes = fes-np.min(fes)
            self.fes = np.array(fes)
            
            
        elif self.cvs == 2:
            cv1min = self.cv1min
            cv2min = self.cv2min
            
            cv1max = self.cv1max
            cv2max = self.cv2max
            
            #if CV is NOT periodic, enlarge the fes
            if not self.periodic[0]:
                cv1range = self.cv1max-self.cv1min
                cv1min -= cv1range*0.15          
                cv1max += cv1range*0.15
            
            if not self.periodic[1]:
                cv2range = self.cv2max-self.cv2min
                cv2min -= cv2range*0.15          
                cv2max += cv2range*0.15
            
            print("Computing CV bin values. ")
            cv1bin = np.ceil((self.cv1-cv1min)*self.res/(cv1max-cv1min))
            cv2bin = np.ceil((self.cv2-cv2min)*self.res/(cv2max-cv2min))
                        
            cv1bin = cv1bin.astype(int)
            cv2bin = cv2bin.astype(int)
            
            s1res = (self.s1[0]*self.res)/(cv1max - cv1min)
            s2res = (self.s2[0]*self.res)/(cv2max - cv2min)
            
            gauss_res = max(8*s1res, 8*s2res)
            gauss_res = int(gauss_res)
            if gauss_res%2 == 0:
                gauss_res += 1
            
            x = np.arange((-(gauss_res-1)/2), ((gauss_res-1)/2)+1)
            gauss = -np.outer(np.exp(-x*x/2.0/s2res/s2res),\
                              np.exp(-x*x/2.0/s1res/s1res))
            
            fes = np.zeros((self.res,self.res))
            for line in range(len(cv1bin)):
                if (line) % 500 == 0:
                    print(f"Constructing free energy surface: {((line+1)/len(cv1bin)):.1%} finished", end="\r")
                
                fes_center = int((self.res-1)/2)
                gauss_center_to_end = int((gauss_res-1)/2)
                #print(f"\ng_res: {gauss_res}, gauss_center_to_end: {gauss_center_to_end}")
                
                fes_to_edit_cv1 = [cv1bin[line]-1-gauss_center_to_end,
                                   cv1bin[line]-1+gauss_center_to_end]
                fes_to_edit_cv2 = [(self.res-cv2bin[line])-gauss_center_to_end,
                                   (self.res-cv2bin[line])+gauss_center_to_end]
                fes_crop_cv1 = [max(0,fes_to_edit_cv1[0]),min(self.res-1,fes_to_edit_cv1[1])]
                fes_crop_cv2 = [max(0,fes_to_edit_cv2[0]),min(self.res-1,fes_to_edit_cv2[1])]
                
                gauss_crop_cv1 = [max(0,gauss_center_to_end-(cv1bin[line]-1)),
                                      gauss_res-1-max(0,(cv1bin[line]-1)+gauss_center_to_end-self.res+1)]
                gauss_crop_cv2 = [max(0,gauss_center_to_end-(self.res-cv2bin[line])),
                                      gauss_res-1-max(0,(self.res-cv2bin[line])+gauss_center_to_end-self.res+1)]
                
                #print(f"\nfes cv1:{fes_crop_cv1}, cv2:{fes_crop_cv2}; \n gauss cv1:{gauss_crop_cv1}, cv2:  {gauss_crop_cv2}")
                
                fes[fes_crop_cv2[0]:fes_crop_cv2[1]+1,fes_crop_cv1[0]:fes_crop_cv1[1]+1]\
                        += gauss[gauss_crop_cv2[0]:gauss_crop_cv2[1]+1,gauss_crop_cv1[0]:gauss_crop_cv1[1]+1]\
                        * self.heights[line]
                
                if self.periodic[0]:
                    if cv1bin[line] < gauss_center_to_end:
                        fes_crop_cv1_p = [self.res-1+(cv1bin[line]-gauss_center_to_end),self.res-1]
                        gauss_crop_cv1_p = [0,gauss_center_to_end-cv1bin[line]]
                        fes[fes_crop_cv2[0]:fes_crop_cv2[1]+1,fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1]\
                                    += gauss[gauss_crop_cv2[0]:gauss_crop_cv2[1]+1,gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1]\
                                    * self.heights[line]
                    
                    if cv1bin[line] > (self.res-gauss_center_to_end):
                        fes_crop_cv1_p = [0,gauss_center_to_end+cv1bin[line]-self.res-1]
                        gauss_crop_cv1_p = [gauss_res-(gauss_center_to_end+cv1bin[line]-self.res),gauss_res-1]
                        fes[fes_crop_cv2[0]:fes_crop_cv2[1]+1,fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1]\
                                    += gauss[gauss_crop_cv2[0]:gauss_crop_cv2[1]+1,gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1]\
                                    * self.heights[line]
                
                if self.periodic[1]:
                    if cv2bin[line] > (self.res-gauss_center_to_end):
                        fes_crop_cv2_p = [self.res+self.res-cv2bin[line]-gauss_center_to_end,self.res-1]
                        gauss_crop_cv2_p = [0,gauss_center_to_end-(self.res-cv2bin[line])-1]
                        fes[fes_crop_cv2_p[0]:fes_crop_cv2_p[1]+1,fes_crop_cv1[0]:fes_crop_cv1[1]+1]\
                                    += gauss[gauss_crop_cv2_p[0]:gauss_crop_cv2_p[1]+1,gauss_crop_cv1[0]:gauss_crop_cv1[1]+1]\
                                    * self.heights[line]
                    
                    if cv2bin[line] <= (gauss_center_to_end):
                        fes_crop_cv2_p = [0,gauss_center_to_end-cv2bin[line]]
                        gauss_crop_cv2_p = [gauss_res-(1+gauss_center_to_end-cv2bin[line]),gauss_res-1]
                        fes[fes_crop_cv2_p[0]:fes_crop_cv2_p[1]+1,fes_crop_cv1[0]:fes_crop_cv1[1]+1]\
                                    += gauss[gauss_crop_cv2_p[0]:gauss_crop_cv2_p[1]+1,gauss_crop_cv1[0]:gauss_crop_cv1[1]+1]\
                                    * self.heights[line]
                
                if self.periodic[0] and self.periodic[1]:
                    if ((cv1bin[line] < gauss_center_to_end) or (cv1bin[line] > (self.res-gauss_center_to_end))) \
                            and ((cv2bin[line] > (self.res-gauss_center_to_end)) or (cv2bin[line] <= (gauss_center_to_end))):
                        fes[fes_crop_cv2_p[0]:fes_crop_cv2_p[1]+1,fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1]\
                                    += gauss[gauss_crop_cv2_p[0]:gauss_crop_cv2_p[1]+1,gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1]\
                                    * self.heights[line]
            print("\n")
            fes = fes-np.min(fes)
            self.fes = np.array(fes)
        elif self.cvs == 3:
            ...
        else:
            print("Fes object doesn't have supported number of CVs.")
    
    def plot(self, png_name=None, contours=True, contours_spacing=0.0, aspect = 1.0, cmap = "jet", 
                 energy_unit="kJ/mol", xlabel=None, ylabel=None, label_size=12, image_size=[10,7], vmin = 0, vmax = None):
        if vmax == None:
            vmax = np.max(self.fes)
            
        if contours_spacing == 0.0:
            default_contours_spacing = (vmax-vmin)/5.0
            contours_spacing = default_contours_spacing
        
        cmap = cm.get_cmap(cmap)
        
        cmap.set_over("white")
        cmap.set_under("white")
        
        if self.cvs >= 1:
            if not self.periodic[0]:
                cv1min = self.cv1min - (self.cv1max-self.cv1min)*0.15
                cv1max = self.cv1max + (self.cv1max-self.cv1min)*0.15
            else:
                cv1min = self.cv1min
                cv1max = self.cv1max 
        if self.cvs >=2:
            if not self.periodic[1]:
                cv2min = self.cv2min - (self.cv2max-self.cv2min)*0.15
                cv2max = self.cv2max + (self.cv2max-self.cv2min)*0.15
            else:
                cv2min = self.cv2min
                cv2max = self.cv2max 
        if self.cvs == 3:
            if not self.periodic[2]:
                cv3min = self.cv3min - (self.cv3max-self.cv3min)*0.15
                cv3max = self.cv3max + (self.cv3max-self.cv3min)*0.15
            else:
                cv3min = self.cv3min
                cv3max = self.cv3max 
        
        if self.cvs == 1:
            plt.figure(figsize=(image_size[0],image_size[1]))
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
            
        if self.cvs == 2:
            plt.figure(figsize=(image_size[0],image_size[1]))
            plt.imshow(self.fes, cmap=cmap, interpolation='nearest', 
                       extent=[cv1min, cv1max, cv2min, cv2max], 
                       aspect = (((cv1max-cv1min)/(cv2max-cv2min))/(aspect)),
                       vmin = vmin, vmax = vmax)
            cbar = plt.colorbar()
            cbar.set_label(energy_unit, size=label_size)
            if contours:
                cont = plt.contour(self.fes, 
                         levels = np.arange(0, (vmax - 0.01), contours_spacing), 
                         extent=[cv1min, cv1max, cv2max, cv2min], 
                         colors = "k")
                plt.clabel(cont, levels = np.arange(0, (vmax - 0.01), contours_spacing))
            if xlabel == None:
                plt.xlabel(f'CV1 - {self.cv1_name}', size=label_size)
            else:
                plt.xlabel(xlabel, size=label_size)
            if ylabel == None:
                plt.ylabel(f'CV2 - {self.cv2_name}', size=label_size)
            else:
                plt.ylabel(ylabel, size=label_size)
        
        if self.cvs == 3:
            ...
        
        
        if png_name != None:
            plt.savefig(name)
        
    def set_fes(self, fes):
        self.fes = fes
        
    def surface_plot(self, cmap = "jet", 
                     energy_unit="kJ/mol", xlabel=None, ylabel=None, zlabel = None, 
                     label_size=12, image_size=[12,7], rstride=1, cstride=1, vmin = 0, vmax = None):
        if self.cvs == 2:
            if not self.periodic[0]:
                cv1min = self.cv1min - (self.cv1max-self.cv1min)*0.15
                cv1max = self.cv1max + (self.cv1max-self.cv1min)*0.15
            else:
                cv1min = self.cv1min
                cv1max = self.cv1max 
            if not self.periodic[1]:
                cv2min = self.cv2min - (self.cv2max-self.cv2min)*0.15
                cv2max = self.cv2max + (self.cv2max-self.cv2min)*0.15
            else:
                cv2min = self.cv2min
                cv2max = self.cv2max
            
            #fes_sp = self.fes
            #if (vmax != None) and vmax < np.max(self.fes):
            #    fes_sp[fes_sp > vmax] = vmax
            # 
            x = np.linspace(self.cv1min, self.cv1max, self.res)
            y = np.linspace(self.cv2max, self.cv2min, self.res)
            
            X, Y = np.meshgrid(x, y)
            Z = self.fes
            
            fig = plt.figure(figsize=(image_size[0],image_size[1]))
            ax = plt.axes(projection='3d')
            ax.plot_surface(X, Y, Z, rstride=rstride, cstride=cstride,
                            cmap=cmap, edgecolor='none', vmin=vmin, vmax=vmax)
            #labels
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
        else:
            print(f"Surface plot only works for FES with exactly two CVs, and this FES has {self.cvs}")
    
    def removeCV(self, CV=None, energy_unit="kJ/mol", temp=300.0):
        CV = int(float(CV))
        print(f"Will remove CV {CV}.")
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
                    new_prob = np.sum(probabilities, axis=1)
                    new_fes = Fes(hills=None)
                    new_fes.fes = -8.314*temp*np.log(new_prob)/1000
                    new_fes.fes = new_fes.fes - np.min(new_fes.fes)
                    new_fes.cvs = 1
                    new_fes.res = self.res
                    new_fes.periodic = [self.periodic[1]]
                    new_fes.cv1min = self.cv2min
                    new_fes.cv1max = self.cv2max
                    new_fes.cv1_name = self.cv2_name
                if CV == 2:
                    new_prob = np.sum(probabilities, axis=0)
                    new_fes = Fes(hills=None)
                    new_fes.fes = -8.314*temp*np.log(new_prob)/1000
                    new_fes.fes = new_fes.fes - np.min(new_fes.fes)
                    new_fes.cvs = 1
                    new_fes.res = self.res
                    new_fes.periodic = [self.periodic[0]]
                    new_fes.cv1min = self.cv1min
                    new_fes.cv1max = self.cv1max
                    new_fes.cv1_name = self.cv1_name
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
                return new_fes
            else:
                print("Error: unknown energy unit")
                return None
        elif self.cvs == 3:
            ...
        
        
    
class Minima(Fes):
    def __init__(self, fes, nbins = 8):
        super().__init__(fes.hills, calculate_new_fes=False)
        self.fes = fes.fes
        self.periodic = fes.periodic
        self.cvs = fes.cvs
        self.res = fes.res

        if self.cvs >= 1:
            self.cv1_name = fes.cv1_name
            self.cv1min = fes.cv1min
            self.cv1max = fes.cv1max
        if self.cvs >= 2:
            self.cv2min = fes.cv2min
            self.cv2max = fes.cv2max
            self.cv2_name = fes.cv2_name
        if self.cvs == 3:
            self.cv3min = fes.cv3min
            self.cv3max = fes.cv3max
            self.cv3_name = fes.cv3_name
        
        self.findminima(nbins=nbins)

    def findminima(self, nbins=8):
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
            if not self.periodic[0]:
                cv1min = self.cv1min - (self.cv1max-self.cv1min)*0.15
                cv1max = self.cv1max + (self.cv1max-self.cv1min)*0.15
            else:
                cv1min = self.cv1min
                cv1max = self.cv1max 
        if self.cvs >=2:
            if not self.periodic[1]:
                cv2min = self.cv2min - (self.cv2max-self.cv2min)*0.15
                cv2max = self.cv2max + (self.cv2max-self.cv2min)*0.15
            else:
                cv2min = self.cv2min
                cv2max = self.cv2max 
        if self.cvs == 3:
            if not self.periodic[2]:
                cv3min = self.cv3min - (self.cv3max-self.cv3min)*0.15
                cv3max = self.cv3max + (self.cv3max-self.cv3min)*0.15
            else:
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
                    min_cv1 = (((min_cv1_b+0.5)/self.res)*(cv1max-cv1min))+cv1min
                    if self.minima == []:
                        self.minima=np.array([round(bin_min, 6), int(min_cv1_b), round(min_cv1, 6)])
                    else:
                        self.minima=np.vstack((self.minima, np.array([round(bin_min, 6), int(min_cv1_b), round(min_cv1, 6)])))
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
                        min_cv1 = (((min_cv1_b+0.5)/self.res)*(cv1max-cv1min))+cv1min
                        if self.minima == []:
                            self.minima=np.array([round(bin_min, 6), int(min_cv1_b), round(min_cv1, 6)])
                        else:
                            self.minima=np.vstack((self.minima, np.array([round(bin_min, 6), int(min_cv1_b), round(min_cv1, 6)])))
            
        elif self.cvs == 2:
            for bin1 in range(0,nbins):
                for bin2 in range(0,nbins):
                    fes_slice = self.fes[bin2*bin_size:(bin2+1)*bin_size,
                                         bin1*bin_size:(bin1+1)*bin_size]
                    bin_min = np.min(fes_slice)
                    argmin = np.argmin(fes_slice)
                    # indexes of global minimum of a bin
                    bin_min_arg_cv1 = int(argmin%bin_size)
                    bin_min_arg_cv2 = int((argmin-bin_min_arg_cv1)/bin_size)
                    # indexes of that minima in the original fes (indexes +1)
                    min_cv1_b = int(bin_min_arg_cv1+bin1*bin_size)
                    min_cv2_b = int(bin_min_arg_cv2+bin2*bin_size)
                    if (bin_min_arg_cv1 > 0 and bin_min_arg_cv1<(bin_size-1)) \
                                    and (bin_min_arg_cv2 > 0 and bin_min_arg_cv2<(bin_size-1)):
                        min_cv1 = (((min_cv1_b+0.5)/self.res)*(cv1max-cv1min))+cv1min
                        min_cv2 = cv2max-(((min_cv2_b+0.5)/self.res)*(cv2max-cv2min))#+cv2min
                        if self.minima == []:
                            self.minima=np.array([round(bin_min, 6), int(min_cv1_b),\
                                                  int(self.res-min_cv2_b), round(min_cv1, 6), round(min_cv2, 6)])
                        else:
                            self.minima=np.vstack((self.minima, np.array([round(bin_min, 6), int(min_cv1_b), \
                                                                          int(self.res-min_cv2_b), round(min_cv1, 6), round(min_cv2, 6)])))
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
                                around.append(self.fes[min_cv2_b_low,min_cv1_b_low])
                            around.append(self.fes[min_cv2_b,min_cv1_b_low])
                            if not(np.isnan(min_cv2_b_high)):
                                around.append(self.fes[min_cv2_b_high,min_cv1_b_low])
                        #1_b
                        if not(np.isnan(min_cv2_b_low)):
                            around.append(self.fes[min_cv2_b_low,min_cv1_b])
                        if not(np.isnan(min_cv2_b_high)):
                            around.append(self.fes[min_cv2_b_high,min_cv1_b])
                        #1_b_high
                        if not(np.isnan(min_cv1_b_high)):
                            if not(np.isnan(min_cv2_b_low)):
                                around.append(self.fes[min_cv2_b_low,min_cv1_b_high])
                            around.append(self.fes[min_cv2_b,min_cv1_b_high])
                            if not(np.isnan(min_cv2_b_high)):
                                around.append(self.fes[min_cv2_b_high,min_cv1_b_high])
                        if bin_min < np.min(around):
                            min_cv1 = (((min_cv1_b+0.5)/self.res)*(cv1max-cv1min))+cv1min
                            min_cv2 = cv2max-(((min_cv2_b+0.5)/self.res)*(cv2max-cv2min))
                            if self.minima == []:
                                self.minima=np.array([round(bin_min, 6), int(min_cv1_b), int(self.res-min_cv2_b), \
                                                      round(min_cv1, 6), round(min_cv2, 6)])
                            else:
                                self.minima=np.vstack((self.minima, np.array([round(bin_min, 6), int(min_cv1_b), \
                                                                              int(self.res-min_cv2_b), round(min_cv1, 6), round(min_cv2, 6)])))
        elif self.cvs == 3:
            ...
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
        elif len(self.minima.shape) == 1:
            self.minima = np.append("A", self.minima)
        
        if self.cvs == 1:
            self.minima = pd.DataFrame(self.minima, columns = ["Minimum", "free energy", "CV1bin", "CV1 - "+self.cv1_name])
        elif self.cvs == 2:
            if len(self.minima.shape)>1:
                self.minima = pd.DataFrame(np.array(self.minima), columns = ["Minimum", "free energy", "CV1bin", "CV2bin", 
                                                               "CV1 - "+self.cv1_name, "CV2 - "+self.cv2_name])
            elif len(self.minima.shape) == 1:
                self.minima = pd.DataFrame([self.minima], columns = ["Minimum", "free energy", "CV1bin", "CV2bin", 
                                                               "CV1 - "+self.cv1_name, "CV2 - "+self.cv2_name])
        elif self.cvs == 3:
            ...
        

    def plot(self, png_name=None, contours=True, contours_spacing=0.0, aspect = 1.0, cmap = "jet", 
                 energy_unit="kJ/mol", xlabel=None, ylabel=None, label_size=12, image_size=[10,7], 
                 color=None, vmin = 0, vmax = None):
        if vmax == None:
            vmax = np.max(self.fes)
        color_set = True
        if color == None:
            color_set = False
        if contours_spacing == 0.0:
            default_contours_spacing = (vmax-vmin)/5.0
            contours_spacing = default_contours_spacing
        
        cmap = cm.get_cmap(cmap)
        
        cmap.set_over("white")
        cmap.set_under("white")
        
        if self.cvs >= 1:
            if not self.periodic[0]:
                cv1min = self.cv1min - (self.cv1max-self.cv1min)*0.15
                cv1max = self.cv1max + (self.cv1max-self.cv1min)*0.15
            else:
                cv1min = self.cv1min
                cv1max = self.cv1max 
        if self.cvs >=2:
            if not self.periodic[1]:
                cv2min = self.cv2min - (self.cv2max-self.cv2min)*0.15
                cv2max = self.cv2max + (self.cv2max-self.cv2min)*0.15
            else:
                cv2min = self.cv2min
                cv2max = self.cv2max 
        if self.cvs == 3:
            if not self.periodic[2]:
                cv3min = self.cv3min - (self.cv3max-self.cv3min)*0.15
                cv3max = self.cv3max + (self.cv3max-self.cv3min)*0.15
            else:
                cv3min = self.cv3min
                cv3max = self.cv3max 
        
        if self.cvs == 1:
            plt.figure(figsize=(image_size[0],image_size[1]))
            X = np.linspace(cv1min, cv1max, self.res)
            plt.plot(X, self.fes)
            
            if not color_set:
                color = "black"
            
            ferange = np.max(self.fes) - np.min(self.fes)
            
            
            if self.minima.shape[0] == 1:
                plt.text(float(self.minima.iloc[0,3]), float(self.minima.iloc[0,1])+ferange*0.05, self.minima.iloc[0,0],
                             fontsize=label_size, horizontalalignment='center',
                             verticalalignment='bottom', c=color)
            elif self.minima.shape[0] > 1:
                for m in range(len(self.minima.iloc[:,0])):
                    plt.text(float(self.minima.iloc[m,3]), float(self.minima.iloc[m,1])+ferange*0.05, self.minima.iloc[m,0],
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
            
        elif self.cvs == 2:        
            plt.figure(figsize=(image_size[0],image_size[1]))
            im = plt.imshow(self.fes, cmap=cmap, interpolation='nearest', 
                       extent=[cv1min, cv1max, cv2min, cv2max], 
                       aspect = (((cv1max-cv1min)/(cv2max-cv2min))/(aspect)),
                       vmin = vmin, vmax = vmax)
            plt.imshow(self.fes, cmap=cmap, interpolation='nearest', 
                       extent=[cv1min, cv1max, cv2min, cv2max], 
                       aspect = (((cv1max-cv1min)/(cv2max-cv2min))/(aspect)),
                       vmin = vmin, vmax = vmax)
            cbar = plt.colorbar()
            cbar.set_label(energy_unit, size=label_size)

            if self.minima.shape[0] == 1:
                background = im.cmap(im.get_array()[int(float(self.minima.iloc[0,2])),int(float(self.minima.iloc[0,3]))]/255)
                luma = background[0]*0.2126+background[1]*0.7152+background[3]*0.0722
                if luma > 0.6 and not color_set:
                    color = "black"
                elif luma <= 0.6 and not color_set:
                    color="white"
                plt.text(float(self.minima.iloc[0,4]), float(self.minima.iloc[0,5]), self.minima.iloc[0,0],
                             fontsize=label_size, horizontalalignment='center',
                             verticalalignment='center', c=color)
            elif self.minima.shape[0] > 1:
                for m in range(len(self.minima.iloc[:,0])):
                    background = im.cmap((float(self.minima.iloc[m,1])-vmin)/(vmax-vmin))
                    luma = background[0]*0.2126+background[1]*0.7152+background[3]*0.0722
                    if luma > 0.6 and not color_set:
                        color = "black"
                    elif luma <= 0.6 and not color_set:
                        color="white"
                    plt.text(float(self.minima.iloc[m,4]), float(self.minima.iloc[m,5]), self.minima.iloc[m,0],
                             fontsize=label_size, horizontalalignment='center',
                             verticalalignment='center', c=color)

            if contours:
                cont = plt.contour(self.fes, 
                         levels = np.arange(0, (vmax+0.01), contours_spacing), 
                         extent=[cv1min, cv1max, cv2max, cv2min], 
                         colors = "k")
                plt.clabel(cont, levels = np.arange(0, (vmax+0.01), contours_spacing))
            if xlabel == None:
                plt.xlabel(f'CV1 - {self.cv1_name}', size=label_size)
            else:
                plt.xlabel(xlabel, size=label_size)
            if ylabel == None:
                plt.ylabel(f'CV2 - {self.cv2_name}', size=label_size)
            else:
                plt.ylabel(ylabel, size=label_size)
        
        elif self.cvs == 3:
            ...
        
        
        if png_name != None:
            plt.savefig(png_name)
        
class FEProfile:
    def __init__(self, minima, hills):
        #self.periodic = minima.periodic
        self.cvs = minima.cvs
        self.res = minima.res
        self.minima = minima.minima
        self.periodic = minima.periodic
        self.heights = hills.get_heights()
        self.sigmas_dont_change = True
        
        if self.cvs >= 1:
            self.cv1_name = minima.cv1_name
            self.cv1min = minima.cv1min
            self.cv1max = minima.cv1max
            self.cv1 = hills.get_cv1()
            self.s1 = hills.get_sigma1()
            if ((np.max(self.s1)/np.min(self.s1))>1.00000001):
                self.sigmas_dont_change = False
        if self.cvs >= 2:
            self.cv2min = minima.cv2min
            self.cv2max = minima.cv2max
            self.cv2_name = minima.cv2_name
            self.cv2 = hills.get_cv2()
            self.s2 = hills.get_sigma2()
            if ((np.max(self.s2)/np.min(self.s2))>1.00000001):
                self.sigmas_dont_change = False
        if self.cvs == 3:
            self.cv3min = minima.cv3min
            self.cv3max = minima.cv3max
            self.cv3_name = minima.cv3_name
            self.cv3 = hills.get_cv3()
            self.s3 = hills.get_sigma3()
            if ((np.max(self.s3)/np.min(self.s3))>1.00000001):
                self.sigmas_dont_change = False
        
        if len(minima.minima.shape)>1:
            self.makefeprofile(hills)
        else: 
            print("There is only one local minimum on the free energy surface.")
        
        
    def makefeprofile(self, hills):
        hillslenght = len(hills.get_cv1())
        
        if hillslenght < 256:
            profilelenght = hillslenght
            scantimes = np.array(range(hillslenght))
        else:
            profilelenght = 256
            scantimes = np.array(((hillslenght/(profilelenght))*np.array((range(1,profilelenght+1)))))
            scantimes -= 1
            scantimes = scantimes.astype(int)
        
        number_of_minima = self.minima.shape[0]
        
        self.feprofile = np.zeros((self.minima.Minimum.shape[0]+1))
               
        if self.sigmas_dont_change:
            if self.cvs == 1:
                cv1min = self.cv1min

                cv1max = self.cv1max

                #if CV is NOT periodic, enlarge the fes
                if not self.periodic[0]:
                    cv1range = self.cv1max-self.cv1min
                    cv1min -= cv1range*0.15          
                    cv1max += cv1range*0.15

                cv1bin = np.ceil((self.cv1-cv1min)*self.res/(cv1max-cv1min))
                cv1bin = cv1bin.astype(int)

                s1res = (self.s1[0]*self.res)/(cv1max - cv1min)

                gauss_res = 8*s1res
                gauss_res = int(gauss_res)
                if gauss_res%2 == 0:
                    gauss_res += 1

                x = np.arange((-(gauss_res-1)/2), ((gauss_res-1)/2)+1)
                gauss = -np.exp(-x*x/2.0/s1res/s1res)

                fes = np.zeros((self.res))

                for line in range(len(cv1bin)):
                    if (line) % 500 == 0:
                        print(f"Constructing free energy surface: {((line+1)/len(cv1bin)):.1%} finished", end="\r")

                    fes_center = int((self.res-1)/2)
                    gauss_center_to_end = int((gauss_res-1)/2)
                    #print(f"\ng_res: {gauss_res}, gauss_center_to_end: {gauss_center_to_end}")

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

                    if (line in scantimes):
                        profileline = [line]
                        for m in range(number_of_minima):
                            profileline.append(fes[int(float(self.minima.iloc[m,2]))]-\
                                               fes[int(float(self.minima.iloc[0,2]))])
                        self.feprofile = np.vstack([self.feprofile, profileline]) 


            elif self.cvs == 2:        
                cv1bin = []
                cv2bin = []

                cv1min = self.cv1min
                cv2min = self.cv2min

                cv1max = self.cv1max
                cv2max = self.cv2max

                #if CV is NOT periodic, enlarge the fes
                if not self.periodic[0]:
                    cv1range = self.cv1max-self.cv1min
                    cv1min -= cv1range*0.15          
                    cv1max += cv1range*0.15

                if not self.periodic[1]:
                    cv2range = self.cv2max-self.cv2min
                    cv2min -= cv2range*0.15          
                    cv2max += cv2range*0.15

                cv1bin = np.ceil((self.cv1-cv1min)*self.res/(cv1max-cv1min))
                cv2bin = np.ceil((self.cv2-cv2min)*self.res/(cv2max-cv2min))

                cv1bin = cv1bin.astype(int)
                cv2bin = cv2bin.astype(int)

                s1res = (self.s1[0]*self.res)/(cv1max - cv1min)
                s2res = (self.s2[0]*self.res)/(cv2max - cv2min)

                gauss_res = max(6*s1res, 6*s2res)
                gauss_res = int(gauss_res)
                if gauss_res%2 == 0:
                    gauss_res += 1

                x = np.arange((-(gauss_res-1)/2), ((gauss_res-1)/2)+1)
                gauss = -np.outer(np.exp(-x*x/2.0/s2res/s2res),\
                                  np.exp(-x*x/2.0/s1res/s1res))
                #print(gauss)
            # crop the unsignificant edges of the gaussian
                #under_threshold_indices = gauss > -1e-4
                #gauss[under_threshold_indices] = 0.0


                fes = np.zeros((self.res,self.res))
                for line in range(len(cv1bin)):
                    if (line) % 500 == 0:
                        print(f"Constructing free energy profile: {((line+1)/len(cv1bin)):.1%} finished", end="\r")

                    fes_center = int((self.res-1)/2)
                    gauss_center_to_end = int((gauss_res-1)/2)
                    #print(f"\ng_res: {gauss_res}, gauss_center_to_end: {gauss_center_to_end}")

                    fes_to_edit_cv1 = [cv1bin[line]-1-gauss_center_to_end,
                                       cv1bin[line]-1+gauss_center_to_end]
                    fes_to_edit_cv2 = [(self.res-cv2bin[line])-gauss_center_to_end,
                                       (self.res-cv2bin[line])+gauss_center_to_end]
                    fes_crop_cv1 = [max(0,fes_to_edit_cv1[0]),min(self.res-1,fes_to_edit_cv1[1])]
                    fes_crop_cv2 = [max(0,fes_to_edit_cv2[0]),min(self.res-1,fes_to_edit_cv2[1])]

                    gauss_crop_cv1 = [max(0,gauss_center_to_end-(cv1bin[line]-1)),
                                          gauss_res-1-max(0,(cv1bin[line]-1)+gauss_center_to_end-self.res+1)]
                    gauss_crop_cv2 = [max(0,gauss_center_to_end-(self.res-cv2bin[line])),
                                          gauss_res-1-max(0,(self.res-cv2bin[line])+gauss_center_to_end-self.res+1)]

                    #print(f"\nfes cv1:{fes_crop_cv1}, cv2:{fes_crop_cv2}; \n gauss cv1:{gauss_crop_cv1}, cv2:  {gauss_crop_cv2}")

                    fes[fes_crop_cv2[0]:fes_crop_cv2[1]+1,fes_crop_cv1[0]:fes_crop_cv1[1]+1]\
                            += gauss[gauss_crop_cv2[0]:gauss_crop_cv2[1]+1,gauss_crop_cv1[0]:gauss_crop_cv1[1]+1]\
                            * self.heights[line]

                    if self.periodic[0]:
                        if cv1bin[line] < gauss_center_to_end:
                            fes_crop_cv1_p = [self.res-1+(cv1bin[line]-gauss_center_to_end),self.res-1]
                            gauss_crop_cv1_p = [0,gauss_center_to_end-cv1bin[line]]
                            fes[fes_crop_cv2[0]:fes_crop_cv2[1]+1,fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1]\
                                        += gauss[gauss_crop_cv2[0]:gauss_crop_cv2[1]+1,gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1]\
                                        * self.heights[line]

                        if cv1bin[line] > (self.res-gauss_center_to_end):
                            fes_crop_cv1_p = [0,gauss_center_to_end+cv1bin[line]-self.res-1]
                            gauss_crop_cv1_p = [gauss_res-(gauss_center_to_end+cv1bin[line]-self.res),gauss_res-1]
                            fes[fes_crop_cv2[0]:fes_crop_cv2[1]+1,fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1]\
                                        += gauss[gauss_crop_cv2[0]:gauss_crop_cv2[1]+1,gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1]\
                                        * self.heights[line]

                    if self.periodic[1]:
                        if cv2bin[line] > (self.res-gauss_center_to_end):
                            fes_crop_cv2_p = [self.res+self.res-cv2bin[line]-gauss_center_to_end,self.res-1]
                            gauss_crop_cv2_p = [0,gauss_center_to_end-(self.res-cv2bin[line])-1]
                            fes[fes_crop_cv2_p[0]:fes_crop_cv2_p[1]+1,fes_crop_cv1[0]:fes_crop_cv1[1]+1]\
                                        += gauss[gauss_crop_cv2_p[0]:gauss_crop_cv2_p[1]+1,gauss_crop_cv1[0]:gauss_crop_cv1[1]+1]\
                                        * self.heights[line]

                        if cv2bin[line] <= (gauss_center_to_end):
                            fes_crop_cv2_p = [0,gauss_center_to_end-cv2bin[line]]
                            gauss_crop_cv2_p = [gauss_res-(1+gauss_center_to_end-cv2bin[line]),gauss_res-1]
                            fes[fes_crop_cv2_p[0]:fes_crop_cv2_p[1]+1,fes_crop_cv1[0]:fes_crop_cv1[1]+1]\
                                        += gauss[gauss_crop_cv2_p[0]:gauss_crop_cv2_p[1]+1,gauss_crop_cv1[0]:gauss_crop_cv1[1]+1]\
                                        * self.heights[line]

                    if self.periodic[0] and self.periodic[1]:
                        if ((cv1bin[line] < gauss_center_to_end) or (cv1bin[line] > (self.res-gauss_center_to_end))) \
                                and ((cv2bin[line] > (self.res-gauss_center_to_end)) or (cv2bin[line] <= (gauss_center_to_end))):
                            fes[fes_crop_cv2_p[0]:fes_crop_cv2_p[1]+1,fes_crop_cv1_p[0]:fes_crop_cv1_p[1]+1]\
                                        += gauss[gauss_crop_cv2_p[0]:gauss_crop_cv2_p[1]+1,gauss_crop_cv1_p[0]:gauss_crop_cv1_p[1]+1]\
                                        * self.heights[line]

                    if (line in scantimes):
                        profileline = [line]
                        for m in range(number_of_minima):
                            profileline.append(fes[(self.res-int(float(self.minima.iloc[m,3]))),int(float(self.minima.iloc[m,2]))]-\
                                               fes[(self.res-int(float(self.minima.iloc[0,3]))),int(float(self.minima.iloc[0,2]))])
                        self.feprofile = np.vstack([self.feprofile, profileline]) 


            elif self.cvs == 3:
                ...
            else:
                print("Fes object doesn't have supported number of CVs.")
        else:
            ...
    
    def plot(self, name="FEprofile.png",image_size=[10,7], xlabel=None, ylabel=None, label_size=12, cmap="jet"):       
        plt.figure(figsize=(image_size[0],image_size[1]))
        
        cmap=cm.get_cmap(cmap)
        
        #colors = cm.jet((self.minima.iloc[:,1].to_numpy()).astype(float)/\
        #                (np.max(self.minima.iloc[:,1].to_numpy().astype(float))))
        colors = cmap(np.linspace(0,1,self.minima.shape[0]))
        for m in range(self.minima.shape[0]):
            plt.plot(self.feprofile[:,0], self.feprofile[:,m+1], color=colors[m])

        if xlabel == None:
            plt.xlabel('time (ps)', size=label_size)
        else:
            plt.xlabel(xlabel, size=label_size)
        if ylabel == None:
            plt.ylabel('free energy difference (kJ/mol)', size=label_size)
        else:
            plt.ylabel(ylabel, size=label_size)
        plt.savefig(name)