from numpy import *
from Plotter import *

class Params():
    def __init__(self):#,small,k):    
        # World
        self.nz = 101 # vertical
        self.nx = 100 # into page
        self.z = linspace(0,1,self.nz)
        if self.nz > 1:
            self.dz = self.z[1] - self.z[0]
        self.Z = tile(self.z,((self.nx,1))).T

        # Time
        self.t = 0.
        self.tstep = 0
        self.marching = 'constant'
        self.t_max = 1.
        self.dt = 0.002
        self.tsteps = int(self.t_max/self.dt)
        self.nslices = 100.
        
        # Initial grainsize distribution
        self.IC = 'bi'
        self.small = 50
        self.s_m = 0.5 # smallest size

        # Loading Parameters
        self.k_s = 1.
        self.k_m = 0.
        self.k_b = 0.
        self.beta = 0.5 # between 0 and 1

        self.shear = 'simple'
        buildColorMap(self)
        # Plotting
        self.superName = 'Segregation_test/dt_' + str(self.dt)
        self.fig_width = 345 # Get this from LaTeX using \showthe\columnwidth
        self.aspect_ratio = (sqrt(5)-1.0)/2.0

    def initialPlots(self,s):
        pass
            
    def finalPlots(self,s,image_store):
        make_s_bar_contour(image_store,self)

    def everyTimestepPlots(self,s):
        pass
    
    def everyNTimestepsPlots(self,s,tstep):
        pass

    def everyNslicesPlots(self,s,tstep):
        pass

    def afterEachShufflePlots(self,s,i,nsims):
        pass
