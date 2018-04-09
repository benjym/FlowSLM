from numpy import *
from Plotter import *

class Params():
    def __init__(self):#,small,k):
        # Plotting
        self.fig_width = 345 # Get this from LaTeX using \showthe\columnwidth
        self.aspect_ratio = (sqrt(5)-1.0)/2.0
    
        # World
        self.nz = 200 # vertical
        self.nx = 200 # into page
        self.z = linspace(0,1,self.nz)
        if self.nz > 1:
            self.dz = self.z[1] - self.z[0]
        self.Z = tile(self.z,((self.nx,1))).T

        # Time
        self.marching = 'constant'
        self.t_max = 5.
        self.dt = 0.0005
        self.tsteps = int(self.t_max/self.dt)
        self.nslices = self.tsteps
        
        # Initial grainsize distribution
        self.IC = 'constant_segregated'
        self.small = 50
        self.s_m = 0.5 # smallest size
        self.seed = 0
        random.seed(self.seed)
        
        # Loading Parameters
        self.k_s = 0.
        self.k_m = 0.05
        self.k_b = 0.
        self.beta = 0.5 # between 0 and 1

        self.shear = 'simple'
        self.superName = 'Diffusion_Example/'#seed_' + str(self.seed)
        buildColorMap(self)

    def initialPlots(self,s):
        pass
            
    def finalPlots(self,s,image_store,d10_store):
#        make_pcolor(image_store,self)
        save_image_store(image_store,self)
        save_d10_store(d10_store,self)

    def everyTimestepPlots(self,s):
        pass
    
    def everyNTimestepsPlots(self,s,tstep):
        pass

    def everyNslicesPlots(self,s,tstep):
        pass

    def afterEachShufflePlots(self,s,i,nsims):
        pass
