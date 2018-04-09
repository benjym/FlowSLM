import sys
from numpy import *
from Plotter import setUpPlots
def printSetupToScreen(P):
    print('k_s = '+str(P.k_s) + ', k_m = '+str(P.k_m)+', k_b = '+str(P.k_b))
    if P.IC == 'fractal':
        print('beta = '+str(P.beta) + ', alpha_i = '+str(P.alpha_i)+', shear = '+P.shear)
    else:
        print('beta = '+str(P.beta) + ', s_m = '+str(P.s_m)+', shear = '+P.shear)
    try:
        print('nz = '+str(P.nz)+', nx = ' + str(P.nx)+', ns = '+str(P.ns))
    except:
        print('nz = '+str(P.nz)+', nx = ' + str(P.nx))
    if P.marching == 'constant':
        min_time = inf
        try:
            try:
                seg_time = nan_to_num(1./P.k_s)*P.dz/(P.s.max()/P.s.min()-1.)
            except:
                seg_time = nan_to_num(1./P.k_s)*P.dz/(1./0.5-1.)
            min_time = minimum(min_time,seg_time)
        except:
            seg_time = inf
        try:
            mix_time = 2./P.k_m*P.dz*P.dz
            min_time = minimum(min_time,mix_time)
        except:
            mix_time = inf
        try:
            break_time = nan_to_num(1./P.k_b)
            min_time = minimum(min_time,break_time)
        except:
            break_time = inf
        print('t_s = '+str(seg_time)+
              ', t_m = ' + str(mix_time)+
              ', t_b = '+str(break_time)+
              ', dt = '+str(P.dt)
              )
        min_time /= 2.
        if P.dt > min_time:
            P.dt = min_time
            P.tsteps = int(P.t_max/P.dt)
            setUpPlots(P)
            print('TIMESTEP PROBLEMS. FORCED dt = ' + str(min_time))
        else:
            print('dt is fine.')
        
def buildWorld(P,O):
    P.t = 0.
    P.tstep = 0
    print('Building world')
    s = zeros((P.nz,P.nx,2)) # size and strength
    if P.IC == 'bi' or P.IC == 'bi_segregated':
        s[:,:,0] = randomiseGrid(P)
    elif P.IC == 'fractal':
        s[:,:,0] = smoothFractal(P)
    elif P.IC == 'constant_segregated':
        s[:,:,0] = smoothGridSorted(P)
    elif P.IC == 'bi_fractal':
        s[:,:,0] = smoothDoubleFractal(P)
    elif P.IC == 'constant':
        s[:,:,0] = randomiseSmoothGrid(P)
#    s[:,:,1] = random.normal(0.5,0.12,((P.nz,P.nx)))
    s[:,:,1] = ones((P.nz,P.nx))
#    s[:,:,1] = random.rand(P.nz,P.nx)
#    print s[:,:,1].max(), s[:,:,1].min()
    printSetupToScreen(P)
    return s
    
def buildShearProfile(P,s_bar=0):
#    if s_bar == 0:
#        s_bar = P.s[0]
    if P.shear == 'simple': # simple shear
        gamma_dot=ones((P.nz,P.nx))
    elif P.shear == 'silbert':  # silbert shear profile vector
        gamma_dot=tile((1-P.z)**0.5,(P.nx,1)).T
    elif P.shear == 'linear':  # linear shear
        gamma_dot=tile((1-P.z),(P.nx,1)).T
    elif P.shear == 'band': # shear band
        gamma_dot=tile(4*(P.z - P.z**2),(P.nx,1)).T
    elif P.shear == 'real':  # real shear vector
        gamma_dot = buildRealShearProfile(P,s_bar)
    elif P.shear == 'none':
        gamma_dot=zeros((P.nz,P.nx))
    try:
        flux_max = P.k_s*gamma_dot.max()*(P.s[-1]/P.s[0]-1.)*P.dt/P.dz # CHECK THIS
        print('Expected max seg probability: ' + str(flux_max))
    except:
        pass
    return P, gamma_dot

def buildRealShearProfile(P,s_bar):
    gamma_dot = sqrt(1.-P.Z)/s_bar
#    if gamma_dot.max()>1. or gamma_dot.min()<0.:
#        print('Error: shear out of bounds. Max is ' +
#            str(gamma_dot.max())+' Min is '+str(gamma_dot.min()))
#        sys.exit(1)
    return gamma_dot

def randomiseGrid(P):
    P.ns = 2
    if P.small == 50:
        P.s = array((P.s_m,P.s_m,P.s_m,
                     P.s_m,P.s_m,1.,1.,1.,1.,1.))
    elif P.small == 40:
        P.s = array((P.s_m,P.s_m,P.s_m,
                     P.s_m,1.,1.,1.,1.,1.,1.))
    elif P.small == 20:
        P.s = array((P.s_m,P.s_m,1.,1.,1.,
                     1.,1.,1.,1.,1.))
    elif P.small == 80:
        P.s = array((P.s_m,P.s_m,P.s_m,P.s_m,
                     P.s_m,P.s_m,P.s_m,P.s_m,1.,1.))
    grid0 = repeat(P.s,P.nx*P.nz/len(P.s)) 
    grid0 = reshape(grid0,(P.nz,P.nx))# - 0.01*random.rand(P.nz,P.nx)
    grid1 = zeros_like(grid0)
    if P.IC == 'bi_segregated':
        return grid0
    if P.IC == 'bi':
        grid = zeros_like(grid0)
        for i in xrange(P.nx): # keep each row? the same
            tt = grid0[:,i]
            random.shuffle(tt)
            grid1[:,i] = tt
        for i in xrange(P.nz): # keep each column? the same
            tt = grid1[i,:]
            random.shuffle(tt)
            grid[i,:] = tt
#         print grid[0,:20]
        return grid

def smoothFractal(P):
    try:
        print('WARNING: used minimum size ' + str(P.s_m))
        F = linspace(0.,1.,P.nx*P.nz)
        d = (F*(1. - P.s_m**(3.-P.alpha_i)) + P.s_m**(3.-P.alpha_i))**(1./(3.-P.alpha_i))
    except:
        F = linspace(1./P.nx,1.,P.nx*P.nz)
        d = F**(1./(3.-P.alpha_i))
    P.s = array([d[0],d[-1]])
    random.shuffle(d)
    d = reshape(d,(P.nz,P.nx))
    return d

def smoothDoubleFractal(P):
    F0 = linspace(0.,1.,(100-P.small)*P.nx*P.nz/100)
    d0 = (F0*(1. - P.s_m**(3.-P.alpha_i)) + P.s_m**(3.-P.alpha_i))**(1./(3.-P.alpha_i))
    F1 = linspace(0.,1.,P.small*P.nx*P.nz/100)
    d1 = P.s_m*(F1*(1. - P.s_m**(3.-P.alpha_i)) + P.s_m**(3.-P.alpha_i))**(1./(3.-P.alpha_i))
    d = concatenate([d0,d1])
    P.s = array([d[0],d[-1]])
    random.shuffle(d)
    d = reshape(d,(P.nz,P.nx))
    return d
    
def randomiseSmoothGrid(P):
    P.ns = P.nx
    P.s = linspace(P.s_m,1.,P.ns)
    grid0 = repeat(P.s,P.nz)
    grid0 = reshape(grid0,(P.nz,P.nx))
    grid = zeros_like(grid0)
    for i in xrange(P.nx): # keep each row? the same
        tt = grid0[:,i]
        random.shuffle(tt)
        grid[:,i] = tt
    return grid
    
def smoothGridSorted(P):
    P.ns = P.nx
    P.s = linspace(P.s_m,1.,P.ns)
    grid0 = repeat(P.s,P.nz)
    grid = reshape(grid0,(P.nz,P.nx))
    return grid
