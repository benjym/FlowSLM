import os
import matplotlib.pyplot as plt
import matplotlib.colors as col
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.colors import LogNorm

from numpy import *

def buildColorMap(P):
    cdict = {'red': ((0.0, 1.0, 1.0),
                     (0.25, 1.0, 1.0),
                     (0.5, 1.0, 1.0),
                     (0.75, 0.902, 0.902),
                     (1.0, 0.0, 0.0)),
             'green': ((0.0, 0.708, 0.708),
                       (0.25, 0.302, 0.302),
                       (0.5, 0.2392, 0.2392),
                       (0.75, 0.1412, 0.1412),
                       (1.0, 0.0, 0.0)),
             'blue': ((0.0, 0.4, 0.4),
                      (0.25, 0.3569, 0.3569),
                      (0.5, 0.6078, 0.6078),
                      (0.75, 1., 1.),
                      (1.0, 1., 1.))}
    P.orange_blue_cmap = col.LinearSegmentedColormap('my_colormap',cdict,256)
    cdict = {'red': ((0.0, 0.0, 0.0),
                     (0.5, 0.5, 0.5),
                     (1.0, 1.0, 1.0)),
             'green': ((0.0, 0.0, 0.0),
                       (0.6, 0.6, 0.5),
                       (1.0, 1.0, 1.0)),
             'blue': ((0.0, 0.0, 0.0),
                      (0.6, 0.6, 0.5),
                      (1.0, 1.0, 1.0))}
    P.my_bone_cmap = col.LinearSegmentedColormap('my_colormap',cdict,256)

def setUpPlots(P):
#    fig_width_pt = 345  
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    fig_width = P.fig_width*inches_per_pt  # width in inches
    fig_height = fig_width*P.aspect_ratio      # height in inches
    fig_size =  [fig_width,fig_height]
    params = {'backend': 'agg',
              'axes.labelsize': 10,
              'text.fontsize': 10,
              'legend.fontsize': 10,
              'xtick.labelsize': 10,
              'ytick.labelsize': 10,
              'text.usetex': True,
              'font.family': 'serif',
              'figure.figsize': fig_size}
    plt.rcParams.update(params)

    if hasattr(P, 'superName'):
#        P.folderName = '/media/long/CA/plots/' + P.superName + '/'
        P.folderName = './plots/' + P.superName + '/'
    else:
        if P.IC == 'bi':
            P.folderName = ('/media/long/CA/plots' +
                          '/breakage_' + str(P.k_b) +
                          '/mix_' + str(P.k_m) +
                          '/seg_'+ str(P.k_s) +
                          '/shear_'+ str(P.shear) + 
                          '/small_'+ str(P.small) +
                          '/beta_'+ str(P.beta) +'/'
                      )
        elif P.IC == 'fractal':
            P.folderName = ('~/media/long/CA/plots' +
                          '/breakage_' + str(P.k_b) +
                          '/mix_' + str(P.k_m) +
                          '/seg_'+ str(P.k_s) +
                          '/shear_'+ str(P.shear) + 
                          '/alpha_'+ str(P.alpha_i) +
                          '/beta_'+ str(P.beta) +'/'
                      )
        elif P.IC == 'constant':
            P.folderName = ('~/media/long/CA/plots' +
                          '/breakage_' + str(P.k_b) +
                          '/mix_' + str(P.k_m) +
                          '/seg_'+ str(P.k_s) +
                          '/shear_'+ str(P.shear) + 
                          '/s_m_'+ str(P.s_m) +
                          '/beta_'+ str(P.beta) +'/'
                      )
    P.folderName = os.path.expanduser(P.folderName)

    if not os.path.isdir(P.folderName):
        os.makedirs(P.folderName)
        print('New folder created')


def make_s_bar_contour(image_store,P):
    #print('Plotting Results')
    plt.figure(1)
    plt.clf()
    if P.nz > 1:
#        print image_store.min(), image_store.max()
        ax = plt.axes([0.1,0.1,0.9,0.9])
        IM = ax.contourf(image_store,
                       origin='lower',
                       extent=[0,P.t_max,0,1],
#                       aspect = P.t_max*(sqrt(5)-1.)/2.,
                       #levels=[1e-4, 1e-3, 1e-2, 1e-1, 1e0],
                       levels = linspace(P.s.min(),1.,10.),
                       #norm=LogNorm(vmin=0.01,vmax=1.),
#                       interpolation=None,
                       #clim=(P.s[0],P.s[-1]),
                       #cmap=plt.cm.jet)
#                       clim=[0.3,1.],
                       cmap=P.orange_blue_cmap)
        CB = plt.colorbar(IM, shrink=.8,
#            ticks=((image_store.min(),image_store.max())))
#            ticks=((P.s[0],P.s[-1])))
            ticks=([0.3,1.]))
        CB.set_label(r'$\bar{s}$',rotation='horizontal',size=10,position=((0,0.8)))
#        ax.set_xticks((0,P.t_max))
#        ax.set_xticklabels((0,P.t_max))
        ax.set_yticks((0,1))
        plt.xlabel('Time, $t$')
        plt.ylabel('Height, $z$')
        fileName = ('IC_' + str(P.IC) +
                    '_nz_' + str(P.nz) +
                    '_nx_' + str(P.nx)
                   )
        plt.savefig(P.folderName + fileName + '.png')#,bbox_inches='tight')
    else:
        plt.plot(linspace(0,P.t_max,P.nslices),image_store.T)
        fileName = ('ns_' +
                    str(P.ns) + '_nz_' + str(P.nz) + '_nx_' + str(P.nx))
        plt.savefig(P.folderName + fileName + '.png')
    print('Wrote files .png at', P.folderName + fileName)

def drawStrengthDistribution(s,P):
    plt.figure()
    plt.clf()
    fileName = 'strength.png'
    plt.hist(s[:,:,1].flatten(),20,normed=True)
    plt.savefig(P.folderName + fileName,dpi=100)
    
def add_temporal_gsd(s,P,fit=False):
    data_lin, edge_lin = histogram(s[:,:,0], bins=100, density=True)
    data_log, edge_log = histogram(s[:,:,0], bins=logspace(-3,0,100), density=True)
    #data_log_n, edge_log_n = histogram(s[:,:,0], bins=logspace(-4.,0.,200))
    x_lin = (edge_lin[:-1] + edge_lin[1:])/2.
    x_log = (edge_log[:-1] + edge_log[1:])/2.
    #data_log_n /= x_log**3.
    #x_log_n = (edge_log_n[:-1] + edge_log_n[1:])/2.
    dx_lin = edge_lin[1:] - edge_lin[:-1]
    dx_log = edge_log[1:] - edge_log[:-1]
    #dx_log_n = edge_log_n[1:] - edge_log_n[:-1]
#    N = (sum(data_log_n) - cumsum(data_log_n))/float(P.nx*P.nz)
    s_sort = sort(s[:,:,0])[0]
    N = cumsum(s_sort[::-1]**-3)[::-1] # number = volume/s**3 \propto s**-3

    fig = plt.figure(99,figsize=[8,6])
    plt.subplot(221)
    plt.plot(x_lin,data_lin)
    plt.ylabel(r'$\phi$')
    plt.xlabel(r'$s$')

    plt.subplot(222)
    plt.plot(x_lin,100.*cumsum(data_lin*dx_lin))
    plt.ylabel(r'\% finer')
    plt.xlabel(r'$s$')

    plt.subplot(223)
    plt.semilogx(x_log,100.*cumsum(data_log*dx_log))
    plt.ylabel(r'\% finer')
    plt.xlabel(r'$s$')

    plt.subplot(224)
    plt.loglog(s_sort,N)
    plt.ylabel(r'$N(\Delta>s)$')
    plt.xlabel(r'$s$')
    #plt.axis([1e-4,1e-1,1e2,1e15])
#    plt.loglog(x_log,100.*cumsum(data_log*dx_log))
#    plt.ylabel(r'\% finer')
#    plt.xlabel(r'$s$')


    if fit:
        if P.IC == 'bi_fractal':
            print('hi')
            min_pt0 = 1.1e-2
            max_pt0 = 3e-2
            min_pt1 = 1.1e-1
            max_pt1 = 3e-1
            pm0 = nonzero(abs(s_sort-min_pt0)==abs(s_sort-min_pt0).min())[0]
            pM0 = nonzero(abs(s_sort-max_pt0)==abs(s_sort-max_pt0).min())[0]
            pm1 = nonzero(abs(s_sort-min_pt1)==abs(s_sort-min_pt1).min())[0]
            pM1 = nonzero(abs(s_sort-max_pt1)==abs(s_sort-max_pt1).min())[0]

            x_fit0 = s_sort[pm0:pM0]
            z0 = polyfit(log(x_fit0),log(N[pm0:pM0]),1)
            N_fit0 = exp(polyval(z0,log(x_fit0)))
            plt.loglog(x_fit0,N_fit0,'k--',linewidth=3)
            F_lin0 = ((x_lin**(3.-z0[0]) - x_lin[0]**(3.-z0[0]))/
             (x_lin[-1]**(3.-z0[0]) - x_lin[0]**(3.-z0[0])))
            F_log0 = ((x_log**(3.-z0[0]) - x_log[0]**(3.-z0[0]))/
             (x_log[-1]**(3.-z0[0]) - x_log[0]**(3.-z0[0])))

            x_fit1 = s_sort[pm1:pM1]
            z1 = polyfit(log(x_fit1),log(N[pm1:pM1]),1)
            N_fit1 = exp(polyval(z1,log(x_fit1)))
            plt.loglog(x_fit1,N_fit1,'k--',linewidth=3)
            F_lin1 = ((x_lin**(3.-z1[0]) - x_lin[0]**(3.-z1[0]))/
             (x_lin[-1]**(3.-z1[0]) - x_lin[0]**(3.-z1[0])))
            F_log1 = ((x_log**(3.-z1[0]) - x_log[0]**(3.-z1[0]))/
             (x_log[-1]**(3.-z1[0]) - x_log[0]**(3.-z1[0])))
            
            F_lin = 0.5*(F_lin0 + F_lin1) # NOT RIGHT!!!
            F_log = 0.5*(F_log0 + F_log1) # NOT RIGHT!!!

            plt.subplot(222)
            plt.plot(x_lin,100.*F_lin,'k--')
            plt.subplot(223)
            plt.semilogx(x_log,100.*F_log,'k--')
            print('alphas = ' + str(-z0[0]) + ', ' + str(-z1[0]))
        else:
            min_pt = 1e-2
            max_pt = 1e-1
            pm = nonzero(abs(s_sort-min_pt)==abs(s_sort-min_pt).min())[0]
            pM = nonzero(abs(s_sort-max_pt)==abs(s_sort-max_pt).min())[0]
            x_fit = s_sort[pm:pM]
            z = polyfit(log(x_fit),log(N[pm:pM]),1)
            N_fit = exp(polyval(z,log(x_fit)))
            plt.loglog(x_fit,N_fit,'k--',linewidth=3)
            P.alpha = -z[0]
            F_lin = ((x_lin**(3.-P.alpha) - x_lin[0]**(3.-P.alpha))/
             (x_lin[-1]**(3.-P.alpha) - x_lin[0]**(3.-P.alpha)))
            F_log = ((x_log**(3.-P.alpha) - x_log[0]**(3.-P.alpha))/
             (x_log[-1]**(3.-P.alpha) - x_log[0]**(3.-P.alpha)))

            plt.subplot(222)
            plt.plot(x_lin,100.*F_lin,'k--')
            plt.subplot(223)
            plt.semilogx(x_log,100.*F_log,'k--')
            print('alpha = ' + str(P.alpha))
    fileName = 'gsd_time.png'
    plt.savefig(P.folderName + fileName,dpi=100)#,bbox_inches='tight')
    
    
def add_coop_plot(s,P,fit=False):
    data_log, edge_log = histogram(s[:,:,0].flatten(), density=True, bins=logspace(-6,0,101))
    x_log = (edge_log[:-1] + edge_log[1:])/2.
    dx_log = edge_log[1:] - edge_log[:-1]

    fig = plt.figure(95,figsize=[8,6])
    ax = fig.add_subplot(111)
    if fit:
        # Shrink current axis by 20%    
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        y_0 = 100.*x_log**(3.-2.)
        y_1 = 100.*x_log**(3.-2.6)
        ax.loglog(x_log,y_0,'g--',label=r'$\alpha=2$')
        ax.loglog(x_log,y_1,'k--',label=r'$\alpha=2.6$')
        plt.title(r'$\beta$ = ' + str(P.beta) + '$, k_m$ = ' + str(P.k_m))
        
        # Put a legend to the right of the current axis
    ax.loglog(x_log,100.*cumsum(data_log*dx_log),label=r'$t=' + str(P.t) + '$')
    ax.set_ylabel(r'\% finer')
    ax.set_xlabel(r'$s$')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # if (P.limit_mode == 'productive_power' or P.limit_mode == 'exp') or (P.limit_mode == 'oded' or P.limit_mode == 'lognormal'):
    # subfolder = P.limit_mode + '/' + P.breakage_mode + '/'
    # fileName = 'n_' + str(P.n).zfill(5) + '_beta_' + str(P.beta).zfill(5) + '_coop.png'
    # if not os.path.isdir('temp/' + subfolder):
    #     os.makedirs('temp/' + subfolder)
    # fileName = subfolder + fileName
    # else:
    plt.ylim([1e-4,1e2])
    plt.xlim([1e-6,1e0])
    fileName = P.limit_mode + '/' + P.breakage_mode + '/' + 'k_m_' + str(P.k_m).zfill(5) + '_coop.png'
    plt.savefig('temp/' + fileName,dpi=100)#,bbox_inches='tight')
    plt.ylim([1,100])
    plt.xlim([1e-3,1])
    fileName = P.limit_mode + '/' + P.breakage_mode + '/' + 'k_m_' + str(P.k_m).zfill(5) + '_coop_zoom.png'
    plt.savefig('temp/' + fileName,dpi=100)#,bbox_inches='tight')
        
def draw_state(s,P,tstep):
    plt.figure(98,figsize=(6,1))
    plt.clf()
    plt.pcolor(s[:,:,0],clim=(0.,1.))
    plt.colorbar()
    fileName = ('tstep_' + str(tstep) + '_shear_' + P.shear + '_ns_' +
                str(P.ns) + '_nz_' + str(P.nz) + '_nx_' + str(P.nx))
    plt.savefig(P.folderName + fileName + '.png', dpi=200)
    
def mesh_gsds(s,P):
    fig = plt.figure()
    plt.clf()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    s = array(s)
    iterations = range(s.shape[0])
    bins = logspace(-2,0,100+1)
    S = zeros((s.shape[0],100))
    for i in range(s.shape[0]):
        s_sort = sort(s[i,:,:,0])[0]
        N = cumsum(s_sort[::-1]**-3)[::-1]
        S[i,:] = N
    X,Y = meshgrid(iterations,bin_centres)
    ax.plot_wireframe(X,Y,S)#, rstride=1, cstride=1,lw=0)
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$s$')
    ax.set_zlabel(r'$N(\Delta>s)/N_{total}$')
    #plt.axis([1e-2,1e0,1e-2,1e0])
    fileName = ('gsd_mesh')
#    plt.show()
    plt.savefig(P.folderName + fileName + '.png', dpi=100)
    
def plot_for_paper(s,P):
    s_sort = sort(s[:,:,0])[0]
    N = cumsum(s_sort[::-1]**-3)[::-1] # number = volume/s**3 \propto s**-3

    fig = plt.figure(10)
    plt.loglog(s_sort,N)
    plt.ylabel(r'$N(\Delta>s)$')
    plt.xlabel(r'$s$')
    #plt.axis([1e-4,1e-1,1e2,1e15])
    fileName = 'for_paper.png'
    plt.savefig(P.folderName + fileName)#,bbox_inches='tight')
    
def add_fitting_line(s,P):
    s_sort = sort(s[:,:,0])[0]
    N = cumsum(s_sort[::-1]**-3)[::-1]
    plt.figure(10)
    min_pt = 1e-4
    max_pt = 1e-1
    pm = nonzero(abs(s_sort-min_pt)==abs(s_sort-min_pt).min())[0]
    pM = nonzero(abs(s_sort-max_pt)==abs(s_sort-max_pt).min())[0]
    x_fit = s_sort[pm:pM]
    z = polyfit(log(x_fit),log(N[pm:pM]),1)
    N_fit = exp(polyval(z,log(x_fit)))
    plt.loglog(x_fit,N_fit,'k--',linewidth=3)
    plt.text(1e-3,1e2,r'$\alpha = '+str(-z[0])+'$',fontsize=6)
    print('alpha = ' + str(-z[0]))
    fileName = 'for_paper_with_fit.png'
    plt.savefig(P.folderName + fileName,bbox_inches='tight')

def plot_number_crushing_events(s,P):
    plt.figure()
    plt.clf()
    plt.hist(s[:,:,1].flatten(),log=True)
    plt.xlabel('Number of crushing events')
    plt.ylabel('Frequency')
    if P.limit_mode == 'productive_power' or P.limit_mode == 'exp':
        subfolder = P.limit_mode + '/' + P.breakage_mode + '/'
        fileName = subfolder + 'n_' + str(P.n).zfill(5) + '_beta_' + str(P.beta).zfill(5) + '_crushing_events.png'
    else:
        fileName = P.limit_mode + '/' + P.breakage_mode + '/' + 'beta_' + str(P.beta).zfill(5) +'_zeta_' + str(P.zeta).zfill(5) + '_k_m_' + str(P.k_m).zfill(5) + '_coop.png'    
    plt.savefig('temp/' + fileName + '.png',bbox_inches='tight')
    # print('Wrote files .png at' + P.folderName + fileName)

def make_pcolor(I,P):
    plt.figure(1)
    plt.clf()
    ax = plt.axes([0.1,0.1,0.9,0.9])
    IM = ax.pcolor(I,
                   #origin='lower',
                   #extent=[0,P.t_max,0,1],
                   #aspect = P.t_max*(sqrt(5)-1.)/2.,
                   #interpolation=None,
                   clim=(P.s[0],P.s[-1]),
                   cmap=P.orange_blue_cmap,
                   edgecolors='w')
#    CB = plt.colorbar(IM, shrink=.8,
##        ticks=((image_store.min(),image_store.max())))
#        ticks=((P.s[0],P.s[-1])))
#    CB.set_label(r'$\bar{s}$',rotation='horizontal',size=10,position=((0,0.8)))
    ax.set_xticks(())
    ax.set_yticks(())
    plt.xlabel('Time, $t$')
    plt.ylabel('Height, $z$')
    fileName = ('ns_' +
                str(P.ns) + '_nz_' + str(P.nz) + '_nx_' + str(P.nx))
    plt.savefig(P.folderName + fileName + '.png',bbox_inches='tight')
    savetxt(P.folderName + fileName + '.txt',I)
    print('Wrote files .png at' + P.folderName + fileName)
    
def save_image_store(I,P):
#    fileName = ('IC_' + str(P.IC) + '_nz_' + str(P.nz) + '_nx_' + str(P.nx))
    fileName=('s_bar')
    savetxt(P.folderName + fileName + '.txt',I)
    print('Wrote file .txt at' + P.folderName + fileName)

def save_d10_store(I,P):
    fileName = ('d10')
    savetxt(P.folderName + fileName + '.txt',I)
    print('Wrote file .txt at' + P.folderName + fileName)
    
def save_s(P,s,time):
    fileName = (str(time) + '_nz_' + str(P.nz) + '_nx_' + str(P.nx) + '_state')
    savetxt(P.folderName + fileName + '.txt',s)
    print('Wrote file .txt at' + P.folderName + fileName)

