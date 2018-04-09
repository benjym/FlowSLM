from numpy import *
from Plotter import *

class Params():
    def __init__(self):#,small,k):
        # World
        self.nz = 101 # vertical
        self.nx = 1000 # into page
        self.z = linspace(0,1,self.nz)
        if self.nz > 1:
            self.dz = self.z[1] - self.z[0]
        self.Z = tile(self.z,((self.nx,1))).T

        # Time
        self.marching = 'constant'
        self.t_max = 100.
        self.t = 0.
        self.tstep = 0
        self.k_m = K_M
        dt_m = self.dz*self.dz/self.k_m/2.
        dt_b = self.dz
        self.dt = minimum(dt_m,dt_b)
        
        self.tsteps = int(self.t_max/self.dt)
        self.nslices = 100.
        self.save_list = array([0.05,0.5,1,2,5,10,25,50,100,250,500,750,1000])
        
        # Initial grainsize distribution
        self.IC = 'fractal'
        self.alpha_i = -2. # 2 = constant over size
        self.seed = SEED
        random.seed(self.seed)
        
        # Loading Parameters
        self.k_s = 0.
        self.k_b = 1.
        self.beta = 0.2 # between 0 and 1
        
        self.shear = 'simple'

        # Plotting
        self.superName = ('breakage_mixing/' +
                          'k_m_' + str(self.k_m) + '/'
                          'seed_' + str(self.seed))
                          
        self.fig_width = 345 # Get this from LaTeX using \showthe\columnwidth
        self.aspect_ratio = (sqrt(5)-1.0)/2.0

    def initialPlots(self,s):
        #plot_for_paper(s,self)
        save_s(self,s[:,:,0],'initial')
    
    def finalPlots(self,s,image_store):
#        add_temporal_gsd(s,P)#,True)
#        plot_for_paper(s,self)
#        add_fitting_line(s,self)
        save_s(self,s[:,:,0],'final')
    
    def everyTimestepPlots(self,s):
        pass
    
    def everyNTimestepsPlots(self,s,tstep):
        if (abs(self.t - self.save_list) < 1e-5).any():
             save_s(self,s[:,:,0],tstep*self.dt)

    def everyNslicesPlots(self,s,tstep):
#            add_temporal_gsd(s,self)
#        draw_state(s,self,tstep)
        pass

    def afterEachShufflePlots(self,s,i,nsims):
        if i == nsims-1:
#            add_temporal_gsd(s,P,True)
            add_fitting_line(s,self)
        elif i%5 == 0:
#            add_temporal_gsd(s,P)
            plot_for_paper(s,self)

