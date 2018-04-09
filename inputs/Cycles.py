from numpy import *
from Plotter import *

class Params():
    def __init__(self):#,small,k):
        # World
        self.nz = 1 # vertical
        self.nx = 100000 # into page
        self.z = linspace(0,1,self.nz)
        if self.nz > 1:
            self.dz = self.z[1] - self.z[0]
        self.Z = tile(self.z,((self.nx,1))).T

        # Time
        self.marching = 'constant'
        self.t = 0.
        self.tstep = 0
        self.t_max = 10
        self.dt = 0.05
        self.tsteps = int(self.t_max/self.dt)
        self.nslices = 100.
        self.save_list = array([1,2,5,10,20,50,100,200,500,1000])
        
        # Initial grainsize distribution
        self.IC = 'fractal'
        self.alpha_i = -2. # 2 = constant over size
        self.seed = 1
        random.seed(self.seed)

        # Loading Parameters
        self.k_s = 0.
        self.k_m = 0.
        self.k_b = 1.
        self.beta = 0.2 # between 0 and 1
#        self.zeta = self.nx # size of local neighbourhood (one-sided)
        self.zeta = 1
        self.shear = 'simple'

        # Plotting
        self.superName = ('breakage/alpha_i_' + str(self.alpha_i) +
                          '/beta_' + str(self.beta) +
                          '/zeta_' + str(self.zeta) +
                          '/seed_' + str(self.seed)
                          )
        self.fig_width = 345 # Get this from LaTeX using \showthe\columnwidth
        self.aspect_ratio = (sqrt(5)-1.0)/2.0

    def initialPlots(self,s):
        save_s(self,s[:,:,0],'initial')
    
    def finalPlots(self,s,image_store):
        pass    

    def everyTimestepPlots(self,s):
        pass
    
    def everyNTimestepsPlots(self,s,tstep):
        pass

    def everyNslicesPlots(self,s,tstep):
        pass

    def afterEachShufflePlots(self,s,i,nsims):
        if (i+1 == self.save_list).any():
            save_s(self,s[:,:,0],str(i+1)+'_final')
