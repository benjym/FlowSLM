#!/usr/bin/python

__doc__ = """
CA_poly.py

This script simulates polydisperse flow down an inclined chute 
for a variety of different shear conditions. Now including: breakage!
"""
__author__ = "benjy"
__version__ = "2.2 31/10/2014"

import sys
import time as timer
from numpy import *
from Plotter import *
from Setup import *
import Operators as O

def timeMarch(s,P,O,gamma_dot):
    image_store=zeros((P.nz,P.nslices))
    d10_store = zeros_like(image_store)
    #print('Beginning time marching')
    tic = timer.clock()
    while P.t <= P.t_max - 1e-5:
        if P.marching == 'constant':
            if P.tstep > P.dt:
                if P.tstep%(int(P.tsteps/P.nslices)) == 0:
                #try:
                    image_store[:,int(P.tstep*P.nslices/P.tsteps)]=sum(s[:,:,0],1)/P.nx
                #except:
                #    pass
                #try:
                    d10_store[:,int(P.tstep*P.nslices/P.tsteps)]=O.get_d10(s,P)
                #except:
                #    pass
                    P.everyNslicesPlots(s,P.tstep)
        elif P.marching == 'dynamic':
            if P.t == 0.:
                s_bar = O.get_s_bar_lr(s,P)
                w = O.make_flux(P,s,gamma_dot,s_bar)
            O.update_timestep(s,P,abs(w).max())
            if P.save_me:
                image_store[:,P.slice]=sum(s[:,:,0],1)/P.nx
                d10_store[:,P.slice]=O.get_d10(s,P)
                P.slice += 1
                if (P.t == P.save_list).any():
                    P.everyNslicesPlots(s,P.t)
        if P.shear == 'real':
            if P.nx == 1:
                s_bar = O.get_s_bar_ud(s)
            else:
                s_bar = O.get_s_bar_lr(s,P)
            gamma_dot = buildRealShearProfile(P,s_bar)
        if P.k_b:
            s = O.breakage(P,s,gamma_dot)
        if P.k_s:
            s_bar = O.get_global_s_bar(s,P)
            w = O.make_flux(P,s,gamma_dot,s_bar)
            for a in range(2):
                s = O.segregate(s,P,a,w)
        if P.k_m:
            for a in range(2):
                s = O.mix(s,P,a)
        P.everyNTimestepsPlots(s,P.tstep)
        if P.t != 0.:
            if P.marching == 'dynamic':
                if P.t%(int(P.t_max/10.)) < 1e-4:
                    toc = timer.clock()
                    print('Done: %1.1f%% Elapsed: %1.4f min. Remaining: %1.4f min.'
                        %(P.t*100./P.t_max, (toc-tic)/60.,
                        ((toc-tic)/P.t*P.t_max-(toc-tic))/60.))
            elif P.marching == 'constant':
                if P.tstep%(int(P.tsteps/10.)) == 0:
                    toc = timer.clock()
                    print('Done: %1.1f%% Elapsed: %1.4f min. Remaining: %1.4f min.'
                        %(P.t*100./P.t_max, (toc-tic)/60.,
                        ((toc-tic)/P.t*P.t_max-(toc-tic))/60.))
        P.tstep += 1
        P.t += P.dt
    return image_store, d10_store, s

################################################################################

def main(P,O,i=0,s=False,small=50):
    #printSetupToScreen(P)
    if s is False:
        s = buildWorld(P,O,small)
    P,gamma_dot = buildShearProfile(P)
    #drawStrengthDistribution(s,P)
    image_store,d10_store,s=timeMarch(s,P,O,gamma_dot)
    P.finalPlots(s,image_store,d10_store)
    print('Done ' + str(i+1) + '!')
    return image_store,d10_store,s

################################################################################

if __name__ == '__main__':
    exec('import inputs.' + sys.argv[1] + ' as ParameterFile')
    try:
        nsims = int(sys.argv[2])
    except:
        nsims = 1
    P = ParameterFile.Params()
    tic = timer.clock()
    setUpPlots(P)
    s = buildWorld(P,O)
    s_store = []
    P.initialPlots(s)
    for i in xrange(nsims):
        image_store,d10_store,s = main(P,O,i,s)
        s_store.append(s)
        random.shuffle(s[0,:,0])
        P.afterEachShufflePlots(s,i,nsims)
#    mesh_gsds(s_store,P)
