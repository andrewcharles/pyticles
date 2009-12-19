#! /usr/local/bin/python

""" 
    Optimised smooth particle implementation.
    Fortran
    Copyright Andrew Charles 2009
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
MAX_STEPS = 10000
NP = 64
XMAX = 8 
YMAX = 8
ZMAX = 8
VMAX = 0.0
dt = 0.05
SPACING = 0.9
LIVE_VIEW = False
SIDE = (4,4,4)
TEMPERATURE = 1.8
HLONG=4.0
HSHORT=2.0

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
    #s.redraw(p)
    #print 'update',time() - t

def redraw(t):
    s.redraw(p)

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
        ,temperature=TEMPERATURE,hlong=HLONG,hshort=HSHORT)

    nl = neighbour_list.SortedVerletList(p,cutoff=4.0)
    #nl_2 = neighbour_list.VerletList(p,cutoff=5.0)
    
    p.nlists.append(nl)
    p.nl_default = nl

    p.forces.append(forces.SpamForce(p,nl))
    p.forces.append(forces.CohesiveSpamForce(p,nl))
    nl.build()
    nl.separations()

    # Use the python spam props to initialise
    spam_properties(p,nl,p.h,p.hlr)

    cnt = 0

def main():
    initialise()
    pyglet.clock.schedule_interval(update,0.05)
    pyglet.clock.schedule_interval(redraw,0.2)
    #pyglet.clock.schedule_interval(s.update_eye,1/2.0)
    pyglet.app.run()

if __name__ == "__main__":
    main()




