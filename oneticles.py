#! /usr/local/bin/python

""" 
    This is a particle app that shows use of the
    particle objects and viewer.

    Copyright Andrew Charles 2008
    All rights reserved.
"""

import sys
import time
import particles
import pyglet
import spyview
import profile
pyglet.options['debug_gl'] = False
from pyglet.gl import *

# Global variables
max_steps = 1000
NP1 = 25 
WINDOW_WIDTH = 640 #640
WINDOW_HEIGHT = 480 #480

p = particles.Particles(NP1)

cnt = 0
fps = 0
tstart = time.time()
rebuild_nl = 1

win = pyglet.window.Window(WINDOW_WIDTH,WINDOW_HEIGHT,caption='Pyticles')
img = pyglet.image.load('p.png')
xmap = WINDOW_WIDTH / float(particles.XMAX)
ymap = WINDOW_HEIGHT / float(particles.YMAX)

def update(dt):
    global cnt,fps,rebuild_nl 
    cnt += 1
    p.update()
    if cnt >= max_steps:
        pyglet.app.exit()

@win.event
def on_draw():
    win.dispatch_events()
    glClear(GL_COLOR_BUFFER_BIT)
    glLoadIdentity()
    win.clear()
    for i in range(p.n):
        img.blit(xmap*p.r[i,0],ymap*p.r[i,1])

@win.event
def on_key_press(symbol,modifiers):
    if symbol == pyglet.window.key.R:
    #    initialise()
        pass
    
def initialise():
    global p,nl_1,nl_2,cnt
    print "Restarting"
    p = particles.Particles(NP1)
    nl_1 = particles.neighbour_list.NeighbourList(p,2.0)
    p.nlists.append(nl_1)
    p.nl_default = nl_1

    nl_1.add_force(particles.forces.HookesForce(p,p,p.nlists[0]))

    for nl in p.nlists:
        nl.build_nl_verlet()
    cnt = 0

def main():
    initialise()
    pyglet.clock.schedule_interval(update,1/20.0)
    pyglet.app.run()
    #loop()

if __name__ == "__main__":
    main()



