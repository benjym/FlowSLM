from numpy import *
from Plotter import *

class Params():
    def __init__(self):#,small,k):
        # Plotting
        self.fig_width = 345 # Get this from LaTeX using \showthe\columnwidth
        self.aspect_ratio = (sqrt(5)-1.0)/2.0
    
        # World
        self.nz = 101 # vertical
        self.nx = 500 # into page
        self.z = linspace(0,1,self.nz)
        if self.nz > 1:
            self.dz = self.z[1] - self.z[0]
        self.Z = tile(self.z,((self.nx,1))).T

        # Time
        self.t = 0
        self.tstep = 0
        self.marching = 'constant'
        self.t_max = 5.
        self.dt = 0.001
        self.tsteps = int(self.t_max/self.dt)
        self.nslices = 100.
        
        # Initial grainsize distribution
        self.IC = 'bi'
        self.small = 50
        self.s_m = 0.5 # smallest size

        # Loading Parameters
        self.k_s = 1.
        self.k_m = 0.05
        self.k_b = 0.
        self.beta = 0.5 # between 0 and 1

        self.shear = 'simple'
        self.superName = ('Diffusion_Seg/k_s_' + str(self.k_s) +
                          '/k_m_' + str(self.k_m) +
                          '/dt_' + str(self.dt))
        buildColorMap(self)

    def initialPlots(self,s):
        pass
            
    def finalPlots(self,s,image_store):
        make_s_bar_contour(image_store,self)
        save_image_store(image_store,self)

    def everyTimestepPlots(self,s):
        pass
    
    def everyNTimestepsPlots(self,s,tstep):
        pass

    def everyNslicesPlots(self,s,tstep):
        pass

    def afterEachShufflePlots(self,s,i,nsims):
        pass
