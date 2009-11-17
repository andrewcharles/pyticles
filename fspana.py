#! /usr/local/bin/python

""" 
    Optimised smooth particle implementation.
    Fortran
    Copyright Andrew Charles 2008
    All rights reserved.
"""

import sys
from time import time
import particles
import c_forces as forces
import pyglet
from pyglet.window import mouse
import pview
import profile
import neighbour_list
from pyglet.gl import *
from properties import spam_properties

# Global variables
MAX_STEPS = 100000
NP = 27
XMAX = 20 
YMAX = 20
ZMAX = 20
VMAX = 0.0
dt = 0.05
SPACING = 2.0
LIVE_VIEW = False
SIDE = (3,3,3)
TEMPERATURE = 0.2

p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit='grid',vmax=VMAX
    ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX)
s = pview.ParticleView(p)

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
    p.update(dt)
    s.redraw(p)
    #print 'update',time() - t
    pass

@s.win.event
def on_draw():
    t = time()
    s.fps = pyglet.clock.get_fps()
    s.clear()
    s.redraw(p)
    #print 'draw',time() - t

@s.win.event
def on_key_press(symbol,modifiers):
    if symbol == pyglet.window.key.R:
        initialise()

def initialise():
    global p,nl_1,nl_2,cnt,buttons
    print "Restarting"
    p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit='grid',vmax=VMAX
        ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX
        ,temperature=TEMPERATURE)

    nl_1 = neighbour_list.VerletList(p,cutoff=2.5)
    nl_2 = neighbour_list.VerletList(p,cutoff=5.0)
    
    p.nlists.append(nl_1)
    p.nlists.append(nl_2)
    p.nl_default = nl_1

    p.forces.append(forces.SpamForce(p,nl_1))
    p.forces.append(forces.CohesiveSpamForce(p,nl_2))

    for nl in p.nlists:
        nl.build()
        nl.separations()

    # Use the python spam props to initialise
    spam_properties(p,nl_1,p.h)

    cnt = 0

def main():
    initialise()
    pyglet.clock.schedule_interval(update,0.05)
    #pyglet.clock.schedule_interval(s.update_eye,1/2.0)
    pyglet.app.run()

if __name__ == "__main__":
    main()




