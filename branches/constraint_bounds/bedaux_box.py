#! /usr/local/bin/python

""" 
    Optimised smooth particle implementation using Fortran
    
    Attempt to use the parameters of Bedaux et al

    Scaling parameters:
    t_s = 562K
    rho_s = 3.12 x 10^4 molm^-3
    x_s = 2.81 nm
    p_s = 200 bar

    At 200 bar, 561.6K
    rho_molar_liq = 3.1 x 10^4 molm^-3
    rho_molar_gas = 0.6
    
    Molar mass of water is 18.0153 gmol-1

    Box size is

    36 x 12 x 12

    of which one 6 x 12 x 12 end is constrained to liquid and
    one 6 x 12 x 12 end is constrained to vapour.

    
    
    Copyright Andrew Charles 2009
    All rights reserved.
"""

SPAMCOMPLETE = True
CFORCES = False
#SPAMCOMPLETE = False
#CFORCES = True

import sys
from time import time
import particles
if SPAMCOMPLETE:
    import spam_complete_force
    import forces
elif CFORCES:
    import c_forces as forces
import pyglet
from pyglet.window import mouse
import pview
import profile
import neighbour_list
from pyglet.gl import *
from properties import spam_properties
import numpy as np
import spkernel, spdensity

# Global variables
MAX_STEPS = 1000

#XMAX = 36 
#YMAX = 12
#ZMAX = 12 
#SIDE = (18,6,6)

XMAX = 9 
YMAX = 3
ZMAX = 3 
SIDE = (9,3,3)
LIQX = 3.0
VAPX = 6.0

VMAX = 0.0
dt = 0.1
SPACING = 1.0
NP = SIDE[0]*SIDE[1]*SIDE[2]
TEMPERATURE = 1.0
HLONG = 4.0
HSHORT = 2.0
RINIT = 'grid'

RHO_LIQ = 1.0
RHO_GAS = 0.2

p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit=RINIT,vmax=VMAX
    ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX)
s = pview.ZPRView(p)
nl = neighbour_list.VerletList(p,cutoff=5.0)
liquid_indices = p.r[:,0] <= LIQX
vapour_indices = p.r[:,0] >= VAPX


lr = p.r[liquid_indices,:].copy()
vr = p.r[vapour_indices,:].copy()

cnt = 0
fps = 0
tstart = time()
rebuild_nl = 1


def set_mass(rho,rho_field,m):
    """ Adjust the masses of a set of sph particles
        to the field density.
        r -- positions
        rho -- densities
    """
    v = m/rho
    m = rho_field * v

def update(t):
    global cnt,fps,rebuild_nl,dt,liquid_indices,vapour_indices,lr,vr
    t = time()
    cnt += 1
    if cnt >= MAX_STEPS:
        pyglet.app.exit()
    else:

        #save the positions of the liquid particles
        #lr = p.r[liquid_indices,:].copy()
        #vr = p.r[vapour_indices,:].copy()
        
        #liquid_indices = p.r[:,0] <= LIQX
        #vapour_indices = p.r[:,0] >= VAPX
        p.update(dt)
        p.r[liquid_indices,:] = lr
        p.r[vapour_indices,:] = vr
        for i in range(liquid_indices.size):
            idx = liquid_indices[i]
            set_mass(p.rho[idx],RHO_LIQ,p.m[idx])
        
        for i in range(vapour_indices.size):
            idx = vapour_indices[i]
            set_mass(p.rho[idx],RHO_GAS,p.m[idx])

    
        #print 'rhomin',p.rho[0:p.n].min()
        #print 'rhomean',p.rho[0:p.n].mean()
        #print 'rhomax',p.rho[0:p.n].max()
    print 'update',time() - t

def redraw(t):
    s.redraw(p)

@s.win.event
def on_draw():
    s.clear()
    s.redraw(p)
    #print 'draw',time() - t

@s.win.event
def on_key_press(symbol,modifiers):
    if symbol == pyglet.window.key.R:
        initialise()

def initialise():
    global p,nl,cnt,buttons,liquid_indices,vapour_indices
    print "Restarting"
    p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit=RINIT,vmax=VMAX
        ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX
        ,temperature=TEMPERATURE,hlong=HLONG,hshort=HSHORT
        ,thermostat_temp=TEMPERATURE,thermostat=True)

    #p.v = (np.random.random([NP,3]) - 0.5) * 2
    liquid_indices = p.r[:,0] <= LIQX
    vapour_indices = p.r[:,0] >= VAPX

    p.m[vapour_indices] = 0.5

    print np.mean(p.t)
    nl = neighbour_list.SortedVerletList(p,cutoff=5.0)
    #nl = neighbour_list.VerletList(p,cutoff=1.0)
    #nl_2 = neighbour_list.VerletList(p,cutoff=5.0)
    
    p.nlists.append(nl)
    p.nl_default = nl

    p.forces.append(spam_complete_force.SpamComplete(p,nl))
    p.forces.append(forces.FortranCollisionForce(p,nl,cutoff=0.6))

    nl.build()
    nl.separations()

    # Use the python spam props to initialise
    #spam_properties(p,nl)
    print 'initial mean temperature',np.mean(p.t)
    print 'initial mean density',np.mean(p.rho)

    cnt = 0

def main():
    initialise()
    pyglet.clock.schedule_interval(update,0.05)
    pyglet.clock.schedule_interval(redraw,0.2)
    #pyglet.clock.schedule_interval(s.update_eye,1/2.0)
    pyglet.app.run()

if __name__ == "__main__":
    main()




