#! /usr/local/bin/python

""" 
    Optimised smooth particle implementation.
    Fortran
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

# Global variables
MAX_STEPS = 100000
XMAX = 20 
YMAX = 20 
ZMAX = 20 
VMAX = 0.0
dt = 0.1
SPACING = 1.0
SIDE = (4,4,4)
NP = SIDE[0]*SIDE[1]*SIDE[2]
TEMPERATURE = 0.2
HLONG = 5.0
HSHORT = 2.5
RINIT = 'grid'

p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit=RINIT,vmax=VMAX
    ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX)
s = pview.ZPRView(p)
nl = neighbour_list.VerletList(p,cutoff=5.0)

cnt = 0
fps = 0
tstart = time()
rebuild_nl = 1

def update(t):
    global cnt,fps,rebuild_nl ,dt
    t = time()
    cnt += 1
    if cnt >= MAX_STEPS:
        pyglet.app.exit()
    else:
        p.update(dt)
        print 'rhomin',p.rho[0:p.n].min()
        print 'rhomean',p.rho[0:p.n].mean()
        print 'rhomax',p.rho[0:p.n].max()
        # Trying to find the cause of instability
        if (p.vdot > 1000000).any():
            print 'isbig'
            pyglet.clock.unschedule(update)
        if np.isnan(p.vdot).any():
            print 'isnan'
            pyglet.clock.unschedule(update)
        if np.isnan(p.v).any():
            print 'isnan'
            pyglet.clock.unschedule(update)
        if (nl.rij[0:nl.nip] < 0.5).any():
            print 'isclose'
            pyglet.clock.unschedule(update)
            print nl.rij[0:nl.nip].min()
    #s.redraw(p)
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
    global p,nl,cnt,buttons
    print "Restarting"
    p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit=RINIT,vmax=VMAX
        ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX
        ,temperature=TEMPERATURE,hlong=HLONG,hshort=HSHORT
        ,thermostat_temp=TEMPERATURE,thermostat=True)

    #p.v = (np.random.random([NP,3]) - 0.5) * 2

    print np.mean(p.t)
    nl = neighbour_list.SortedVerletList(p,cutoff=5.0)
    #nl = neighbour_list.VerletList(p,cutoff=1.0)
    #nl_2 = neighbour_list.VerletList(p,cutoff=5.0)
    
    p.nlists.append(nl)
    p.nl_default = nl

    if SPAMCOMPLETE:
        p.forces.append(spam_complete_force.SpamComplete(p,nl))
        #p.forces.append(forces.CollisionForce3d(p,nl,cutoff=0.6))
        p.forces.append(forces.FortranCollisionForce(p,nl,cutoff=0.6))
    if CFORCES:
        p.forces.append(forces.SpamForce(p,nl))
        p.forces.append(forces.CohesiveSpamForce(p,nl))
        p.forces.append(forces.SpamConduction(p,nl))

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




