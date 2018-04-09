import sys
from numpy import *

def update_timestep(s,P,w):
    if P.marching == 'constant':
        pass
    elif P.marching == 'dynamic':
        if P.k_m:
            mix_time = P.dz*P.dz/P.k_m
        else:
            mix_time = 1000
        if P.k_s:
            try:
                w = minimum(P.set_sat_lim,w) # let things go above limit!
            except:
                pass
            seg_time = P.dz/w
        else:
            seg_time = 1000
        if P.k_b:
            break_time = 1./P.k_b
        else:
            break_time = 1000
        P.dt = minimum(minimum(mix_time,seg_time),break_time)
        P.dt = minimum(10.,P.dt) # limit to dt=10.
        P.save_me = 0
        if P.t == 0.: # save first step
            P.save_me = 1
        if (P.t == P.t_list).any():
            P.save_me = 1 
        else:
            skipped_list = (P.t_list>P.t)*(P.t_list<P.dt+P.t)
            if skipped_list.any():
                next_time = nonzero(skipped_list)[0][0]
#                print('Moving to ' + str(P.t_list[next_time]))
                P.dt = P.t_list[next_time] - P.t
#                if trial > 0 and trial < 1.:
#                P.dt = trial
#                P.save_me = 1
#        print break_time, seg_time, mix_time, P.dt, P.t

def swap(s,P,z,dz,x):
    if z == 0:
        dz = maximum(0,dz)
    elif z == P.nz-1:
        dz = minimum(0,dz)
    s[z,x], s[z+dz,x] = s[z+dz,x].copy(), s[z,x].copy()
    return s

def swap_if_size(s,P,z,dz,x):
    if z == 0:
        dz = maximum(0,dz)
    elif z == P.nz-1:
        dz = minimum(0,dz)
    if dz > 0 and s[z,x,0] > s[z+dz,x,0]:
        s[z,x], s[z+dz,x] = s[z+dz,x].copy(), s[z,x].copy()
    elif dz < 0 and s[z,x,0] < s[z+dz,x,0]:
        s[z,x], s[z+dz,x] = s[z+dz,x].copy(), s[z,x].copy()
    return s


def make_flux(P,s,gamma_dot,s_bar):
    w = P.k_s*gamma_dot*(s[:,:,0]/s_bar-1.)
    if P.marching == 'constant':
        flux = w*P.dt/P.dz
        if abs(flux).max()>1.:
            print('flux is out of range. Max is: ' + str(amax(flux)) +
                  ' Min is: ' + str(amin(flux)))
            print('k is '+str(P.k_s)+' gamma_dot is '+str(amax(gamma_dot)))
    return w

def get_d10(s,P):
    d_10 = zeros([P.nz,])
    if P.k_b > 0.:
        bins=logspace(-8,0,101)
    else:
        bins=P.nx/5
    for i in range(P.nz):
        phi, edges = histogram(s[i,:],
                               bins=bins,
                               density=True)

        mid = (edges[:-1] + edges[1:])/2.
        dx = edges[1:] - edges[:-1]
        F = cumsum(phi*dx)
        id_10 = (abs(F-0.1)).argmin()
        d_10[i] = mid[id_10]
#        print d_10
    return d_10

def get_global_s_bar(s,P):
    s_bar = tile(sum(s[:,:,0],1)/P.nx,(P.nx,1)).T
    return s_bar
    
def get_s_bar_lr(s,P):
    # try:
    if P.zeta_mode == "constant":
        if P.zeta == 1:
            s_bar = column_stack(((s[:,-1,0] + s[:,1,0])/2.,
                                  (s[:,:-2,0] + s[:,2:,0])/2.,
                                  (s[:,-2,0] + s[:,0,0])/2.))
        elif P.zeta == P.nx/2:
            s_bar = tile(sum(s[:,:,0],1)/P.nx,(P.nx,1)).T
        elif P.zeta <= 500:
            s_bar = zeros_like(s[:,:,0])
            for i in range(1,P.zeta+1):
                temp1 = roll(s[:,:,0], i,axis=1)
                temp2 = roll(s[:,:,0],-i,axis=1)
                s_bar += temp1 + temp2
            s_bar /= P.zeta*2.
    elif P.zeta_mode == 'variable':
        s_bar = zeros_like(s[:,:,0])
        s_glob = sum(s[:,:,0],1)/P.nx
        # zeta_0 = 100.
        for i in xrange(P.nz):
            for j in xrange(P.nx):
                # zeta = int(around(P.nx*s[i,j,0]/2.))
                # zeta = int(around((P.zeta-1)*s[i,j,0]**3. + 1.))
                zeta = int(ceil(1.5*(s[i,j,0]/s_glob[i])**2))
                # print zeta
                # zeta = int(around((P.zeta-1.)*s[i,j,0] + 1.))
                s_bar[i,j] = get_local_s_bar_lr(s,P,i,j,zeta)
    elif P.zeta_mode == "by_number":
        s_bar = column_stack(((s[:,-1,0]**-2 + s[:,1,0]**-2)/(s[:,-1,0]**-3 + s[:,1,0]**-3),
                              (s[:,:-2,0]**-2 + s[:,2:,0]**-2)/(s[:,:-2,0]**-3 + s[:,2:,0]**-3),
                              (s[:,-2,0]**-2 + s[:,0,0]**-2)/(s[:,-2,0]**-3 + s[:,0,0]**-3)))
    else:
        s_bar = zeros_like(s[:,:,0])
        for i in xrange(P.nx):
            for j in xrange(1,P.zeta+1):
                s_bar[:,i] += s[:,i+j,0] + s[:,i-j,0]
    # except: # if P.zeta not defined
    #     s_bar = column_stack(((s[:,-1,0] + s[:,1,0])/2.,
    #                           (s[:,:-2,0] + s[:,2:,0])/2.,
    #                           (s[:,-2,0] + s[:,0,0])/2.))
    return s_bar

def get_local_s_bar_lr(s,P,i,j,zeta):
    # print s[i,j,0], zeta
    if j+zeta > P.nx-1:
        s_bar = mean(hstack([s[i,j-zeta:j,0], s[i,j+1:,0], s[i,:j+zeta-P.nx+1,0]]))
    elif j-zeta < 0:
        s_bar = mean(hstack([s[i,j-zeta:,0], s[i,:j,0], s[i,j+1:j+zeta+1,0]]))
    else:
        s_bar = (mean(s[i,j-zeta:j,0]) + mean(s[i,j+1:j+1+zeta,0]))/2.
    return s_bar
    
def get_neighbours(P,j,zeta):
    if j+zeta > P.nx-1:
        n = range(j-zeta,j) + range(j+1,P.nx) + range(0,j+zeta-P.nx+1)
        if len(n) != 2*zeta:
            print P.nx, j, zeta
            print 'a', n, len(n), len(range(j-zeta,j)), len(range(j+1,P.nx)), len(range(0,j+zeta-P.nx+1))
            sys.exit()
    elif j-zeta < 0:
        n = range(P.nx+j-zeta,P.nx) + range(0,j) + range(j+1,j+zeta+1)
        if len(n) != 2*zeta:
            print P.nx, j, zeta
            print 'b', n, len(n), len(range(P.nx-j-zeta,P.nx)), len(range(0,j)), len(range(j+1,j+zeta+1))
            sys.exit()
    else:
        n = range(j-zeta,j) + range(j+1,j+1+zeta)
    # print j, n
    return n

def get_s_bar_lr_below(s):
    s_bar = column_stack(((s[:,-1,0] + s[:,1,0])/2.,
                          (s[:,:-2,0] + s[:,2:,0])/2.,
                          (s[:,-2,0] + s[:,0,0])/2.))
    s_bar = row_stack((s_bar[-1,:],s_bar[0:-1,:]))
    return s_bar


def get_s_bar_ud(s):
    #s_bar = tile(sum(s[:,:,0],1)/P.nx,(P.nx,1)).T
    s_bar = row_stack(((s[-1,:,0] + s[1,:,0])/2.,
                          (s[:-2,:,0] + s[2:,:,0])/2.,
                          (s[-2,:,0] + s[0,:,0])/2.))
    return s_bar

def segregate(s,P,a,w):
    for x in range(P.nx):
        for z in range(a,P.nz,2):
            if random.rand(1)<abs(w[z,x]*P.dt/P.dz):
                u_or_d = int(sign(w[z,x]))
                swap_if_size(s,P,z,u_or_d,x)
    return s

def mix(s,P,a):
    u_or_d = sign(random.rand(P.nz,P.nx)-0.5)
    for x in range(P.nx):
        for z in range(a,P.nz,2):
            if random.random()<P.k_m*P.dt/P.dz/P.dz:
                swap(s,P,z,u_or_d[z,x],x)
#                swap(s,P,z,-1,x) # ONLY DOWN
    return s

def breakage(P,s,gamma_dot):
    s_bar = get_s_bar_lr(s,P)
    for z in range(P.nz):
        for x in range(P.nx):
            if P.limit_mode == 'original':
                if abs(s_bar[z,x]-s[z,x,0]) < P.beta*s[z,x,0]:#*s[z,x,1]:
                    s = actual_breakage(P,s,gamma_dot,z,x)
            elif P.limit_mode == 'original_any':
                n = get_neighbours(P,x,P.zeta)
                if (abs(s[z,n,0]-s[z,x,0]) < P.beta*s[z,x,0]).any():#*s[z,x,1]:
                    s = actual_breakage(P,s,gamma_dot,z,x)
            elif P.limit_mode == 'productive':
                if s[z,x,0]/s_bar[z,x] > (1.+P.beta) and s[z,x,0]/s_bar[z,x] < 1./(1.+P.beta):
                    s = actual_breakage(P,s,gamma_dot,z,x)
            elif P.limit_mode == 'productive_any':
                n = get_neighbours(P,x,P.zeta)
                if ((s[z,x,0]/s[z,n,0] < 1.+P.beta)*(s[z,x,0]/s[z,n,0] > 1./(1.+P.beta))).any():
                    s = actual_breakage(P,s,gamma_dot,z,x)
            elif P.limit_mode == 'productive_power': # beta between 0 and 1
                if (3.*abs(s[z,x,0]/s_bar[z,x])**2 - 1.)**P.n < P.beta*(s[z,x,0]**(3./6.)):
                    s = actual_breakage(P,s,gamma_dot,z,x)
            elif P.limit_mode == 'exp': 
                if exp((s[z,x,0] - s_bar[z,x])**2/(2.*(P.n**2.))) < P.beta*(s[z,x,0]**(3./6.)):
                    s = actual_breakage(P,s,gamma_dot,z,x)
            elif P.limit_mode == 'oded':
                if (s[z,x,0]/s_bar[z,x])**2*exp((s[z,x,0] - s_bar[z,x])**2/(2.*(P.n**2.))) < P.beta*(s[z,x,0]**(3./6.)):
                    s = actual_breakage(P,s,gamma_dot,z,x)
            elif P.limit_mode == 'lognormal': 
                if exp(log(s[z,x,0]/s_bar[z,x])**2/(2.*(P.n**2.))) < P.beta*(s[z,x,0]**(3./6.)):
                    s = actual_breakage(P,s,gamma_dot,z,x)
                    # n IS APPROX NUMBER OF DECADES OF SPREAD!
            elif P.limit_mode == 'simple':
                if abs(s_bar[z,x]-s[z,x,0]) < P.beta:
                    s = actual_breakage(P,s,gamma_dot,z,x)
            else: sys.exit('Unknown limit mode:' + str(P.limit_mode))
    return s    

def actual_breakage(P,s,gamma_dot,z,x):
    if random.random() < P.k_b*gamma_dot[z,x]*P.dt:
        if P.breakage_mode == 'normal':
            s[z,x,0] *= random.random() # 'normal'
        elif P.breakage_mode == 'limited':
            s[z,x,0] *= 0.0001/s[z,x,0] + (1.-0.0001/s[z,x,0])*(random.random()**1.) # minimum size limited
        elif P.breakage_mode == 'cascade':e0
            s[z,x,0] *= 0.5 + 0.5*random.random() # cascade
        elif P.breakage_mode == 'lognormal':
            s[z,x,0] *= minimum(1.,random.lognormal(-0.75,0.2))
        elif P.breakage_mode == 'constant':
            s[z,x,0] *= 0.5
        else:
            sys.quit('No breakage mode.')
        s[z,x,1] += 1
    return s
    