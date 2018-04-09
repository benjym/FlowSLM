from numpy import *
from Plotter import *

class Params():
    def __init__(self):#,small,k):
        # World
        self.nz = 10 # vertical
        self.nx = 1000 # into page
        self.z = linspace(0,1,self.nz)
        if self.nz > 1:
            self.dz = self.z[1] - self.z[0]
        self.Z = tile(self.z,((self.nx,1))).T

        # Time
        self.marching = 'constant'
        self.t = 0.
        self.tstep = 0
        self.t_max = 50.
        self.dt = 0.05
        self.tsteps = int(self.t_max/self.dt)
        self.nslices = 100.
        self.save_list = array([0.05,0.5,5.])
        
        # Initial grainsize distribution
        self.IC = 'bi_fractal'
        self.alpha_i = -2.
        self.s_m = MINSIZE # 0.1
        self.small = SMALL # 50
        random.seed = 1

        # Loading Parameters
        self.k_s = 0.
        self.k_m = MIXING
        self.k_b = 1.
        self.beta = 0.2 # between 0 and 1
        
        self.alpha = 2.
        self.shear = 'simple'

        # Plotting
        self.superName = ('breakage/s_m_' + str(self.s_m) +
                          '/small_' + str(self.small) + 
                          '/beta_' + str(self.beta) + 
                          '/D_' + str(self.k_m))
        self.fig_width = 345 # Get this from LaTeX using \showthe\columnwidth
        self.aspect_ratio = (sqrt(5)-1.0)/2.0

    def initialPlots(self,s):
#         plot_for_paper(s,self)
        add_temporal_gsd(s,self)
#         save_s(self,s[:,:,0],'initial')
    
    def finalPlots(self,s,image_store):
#        add_temporal_gsd(s,P)#,True)
#        plot_for_paper(s,self)
#        add_fitting_line(s,self)
        save_s(self,s[:,:,0],'final')
        add_temporal_gsd(s,self)

    
    def everyTimestepPlots(self,s):
        pass
    
    def everyNTimestepsPlots(self,s,tstep):
        if (tstep*self.dt == self.save_list).any():
             save_s(self,s[:,:,0],tstep*self.dt)
             add_temporal_gsd(s,self)

    def everyNslicesPlots(self,s,tstep):
#         add_temporal_gsd(s,self)
#         draw_state(s,self,tstep)
        pass

    def afterEachShufflePlots(self,s,i,nsims):
#         if i == nsims-1:
# #            add_temporal_gsd(s,P,True)
#             add_fitting_line(s,self)
#         elif i%5 == 0:
# #            add_temporal_gsd(s,P)
#             plot_for_paper(s,self)
        pass
