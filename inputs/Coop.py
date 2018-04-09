from numpy import *
from Plotter import *

class Params():
    def __init__(self):#,small,k):
        # World
        self.nz = 21 # vertical
        self.nx = 1000 # into page
        self.z = linspace(0,1,self.nz)
        if self.nz > 1:
            self.dz = self.z[1] - self.z[0]
        self.Z = tile(self.z,((self.nx,1))).T

        # Time
        self.marching = 'constant'
        self.t_max = 1000.
        self.t = 0.
        self.tstep = 0
        self.k_m = K_M
        self.k_s = 0.
        self.k_b = 1.
        
        dt_b = 0.1*self.k_b
        try:
            dt_m = self.dz*self.dz/self.k_m/2.
            self.dt = minimum(dt_m,dt_b)
        except:
            self.dt = dt_b
            dt_m = 0.

        
        self.tsteps = int(self.t_max/self.dt)
        self.nslices = 100.
        self.save_list = array([0.1,0.3,1,3,10,30,100,300,1000,3000,10000])
        # self.save_list = around(logspace(0,5,11))/10. # to nearest 0.1
        # print self.save_list
        # Initial grainsize distribution
        self.IC = 'fractal'
        self.alpha_i = 0. # 2 = constant over size
        # self.IC = 'constant'
        # self.s_m = 0.7
        self.seed = 0
        random.seed(self.seed)
        
        # Loading Parameters
        self.n = 0.2
        self.beta = 10**6
        self.breakage_mode = "BREAKAGE_MODE" #'lognormal' # 'normal'
        self.limit_mode = "lognormal"
        self.shear = 'simple'
        self.zeta_mode = "by_number" #"constant" # "variable"
        self.zeta = 1.

        # Plotting
        self.superName = ('coop/' + 
                          'limit_mode_' + self.limit_mode + '/'
                          'breakage_mode_' + self.breakage_mode + '/' +
                          'beta_' + str(self.beta) + '/' +
                          'alpha_i_' + str(self.alpha_i) + '/' +
                          'k_m_' + str(self.k_m))
                          
        self.fig_width = 345 # Get this from LaTeX using \showthe\columnwidth
        self.aspect_ratio = (sqrt(5)-1.0)/2.0

    def initialPlots(self,s):
        add_coop_plot(s,self,True)
        # save_s(self,s[:,:,0],'initial')
    
    def finalPlots(self,s,image_store,d10_store):
        add_coop_plot(s,self)
        # plot_number_crushing_events(s,self)
        # save_s(self,s[:,:,0],'final')
        # pass
        
    def everyTimestepPlots(self,s):
        pass
    
    def everyNTimestepsPlots(self,s,tstep):
        if (abs(self.t - self.save_list) < 1e-5).any():
             # save_s(self,s[:,:,0],tstep*self.dt)
             add_coop_plot(s,self)
             

    def everyNslicesPlots(self,s,tstep):
        pass

    def afterEachShufflePlots(self,s,i,nsims):
        pass

