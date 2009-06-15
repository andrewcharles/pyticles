#! /usr/local/bin/python

""" 
    The is the front end for the python SPH code most similar
    to what I am doing in the fortran code.

    This is the test bed for new algorithms.

    Copyright Andrew Charles 2008
    All rights reserved.
"""

import sys
import time
import particles
import c_forces as forces
import pyglet
from pyglet.window import mouse
import pview
import profile
import neighbour_list
from pyglet.gl import *

# Global variables
MAX_STEPS = 1000
NP1 = 10

p = particles.SmoothParticleSystem(NP1,maxn=NP1,d=3)
s = pview.ParticleView()

cnt = 0
fps = 0
tstart = time.time()
rebuild_nl = 1

def update(dt):
    global cnt,fps,rebuild_nl 
    cnt += 1
    print 'Step',cnt
    if cnt >= MAX_STEPS:
        pyglet.app.exit()
    p.update()
    pass

@s.win.event
def on_draw():
    s.fps = pyglet.clock.get_fps()
    s.clear()
    s.redraw(p)

@s.win.event
def on_key_press(symbol,modifiers):
    if symbol == pyglet.window.key.R:
        initialise()

def initialise():
    global p,nl_1,nl_2,cnt,buttons
    print "Restarting"
    p = particles.SmoothParticleSystem(NP1,maxn=NP1,d=3)
    
    nl_1 = neighbour_list.VerletList(p,cutoff=2.0)
    nl_2 = neighbour_list.VerletList(p,cutoff=5.0)
    
    p.nlists.append(nl_1)
    p.nlists.append(nl_2)
    p.nl_default = nl_1

    p.forces.append(forces.SpamForce(p,nl_1))
    p.forces.append(forces.CohesiveSpamForce(p,nl_2))

    for nl in p.nlists:
        nl.build()
    cnt = 0

def main():
    initialise()
    pyglet.clock.schedule_interval(update,1/5.0)
    pyglet.app.run()

if __name__ == "__main__":
    main()




