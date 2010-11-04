#! /usr/local/bin/python

""" 
    nanobox_quench

    Optimised smooth particle implementation using Fortran

    Starting with one of the equilibrated states, run a quench.
    
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

    The initial plan was to have one 6 x 12 x 12 end 
    constrained to liquid and one 6 x 12 x 12 end 
    constrained to vapour. I am no longer doing this.

    Copyright Andrew Charles 2009
    All rights reserved.

    You know, this script should define its own viewer...
    That would stop the particle viewer becoming a cruft of
    options.


"""

SPAMCOMPLETE = True

import sys
from time import time
import particles
import spam_complete_force
import forces
import pyglet
import box
from pyglet.window import mouse
import pview
import neighbour_list
from pyglet.gl import *
from properties import spam_properties
import numpy as np
import spkernel, spdensity

# Global variables
MAX_STEPS = 100000

XMAX = 8
YMAX = 8
ZMAX = 8 
NDIM = 3
#SIDE = (10,10,10)
SIDE = (7,7,7)
VMAX = 0.0
dt = 0.01
SPACING = 0.5
NP = SIDE[0]*SIDE[1]*SIDE[2]
TEMPERATURE = 0.8
HLONG = 3.0
HSHORT = 1.5
RINIT = 'grid'
ascl = 7.45e+04
bscl = 5.84e-01
kbscl = 3.29e+04
pmass = 1.386e-01

cnt = 0
fps = 0

box = box.PeriodicBox(xmax=XMAX,ymax=YMAX,zmax=ZMAX)
#box = box.MirrorBox(xmax=XMAX,ymax=YMAX,zmax=ZMAX)
p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit=RINIT,vmax=VMAX
    ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX
    ,temperature=TEMPERATURE,hlong=HLONG,hshort=HSHORT,
    thermostat_temp=TEMPERATURE,thermostat=True,mass=pmass,simbox=box)
nl = neighbour_list.VerletList(p,cutoff=4.0)
s = pview.ZPRView(p)
p.nlists.append(nl)
p.nl_default = nl
p.forces.append(spam_complete_force.SpamComplete(p,nl,adash=ascl,bdash=bscl
    ,kbdash=kbscl,cgrad=0.0,eta=0.0,zeta=0.0))
nl.build()
nl.separations()
spam_properties(p,nl)

print 'initial mean temperature',np.mean(p.t)
print 'initial mean density',np.mean(p.rho)

cnt = 0
fps = 0
tstart = time()
rebuild_nl = 1

# CHANGE THIS CODE TO REDRAW ONLY ON UPDATE
#

def update(t):
    global cnt,fps,rebuild_nl,dt,liquid_indices,vapour_indices,lr,vr
    t = time()
    cnt += 1
    if cnt >= MAX_STEPS:
        pyglet.app.exit()
    else:
        p.update(dt)
        #print 'update',time() - t
        #if np.any(nl.drij[:,:] > 2.0):
        #    print 'wtf'
        #    print (nl.drij > 2.0).size

def redraw(t):
    s.redraw(p)

@s.win.event
def on_draw():
    t = time()
    s.clear()
    s.redraw(p)
    #print 'draw',time() - t


def main():
    pyglet.clock.schedule_interval(update,0.05)
    pyglet.clock.schedule_interval(redraw,0.2)
    #pyglet.clock.schedule_interval(s.update_eye,1/2.0)
    pyglet.app.run()

if __name__ == "__main__":
    main()




