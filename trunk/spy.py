#! /usr/local/bin/python

""" This is the main file for the python sph code
    Write for 3 dimensions as much as possible
    Copyright Andrew Charles 2008
    All rights reserved.
"""

import sys
#sys.path.append("/home/ac/pyglet-1.0")
import time
import particles
import pyglet
import spyview
import profile

# Global variables
max_steps = 1000

p = particles.Particles()
nl_1 = particles.neighbour_list.NeighbourList(p,3.0)
nl_2 = particles.neighbour_list.NeighbourList(p,6.0)

s = spyview.ParticleView()

cnt = 0
fps = 0
tstart = time.time()
rebuild_nl = 1

def initialise():
    global p,nl_1,nl_2,cnt
    print "Restarting"
    p = particles.Particles()
    # a long and a short neighbour list
    nl_1 = particles.neighbour_list.NeighbourList(p,5.0)
    nl_2 = particles.neighbour_list.NeighbourList(p,10.0)
    p.nlists.append(nl_1)
    p.nlists.append(nl_2)
    p.nl_default = nl_1

    p.add_force(particles.forces.SpamForce(p,p.nlists[0]))
    p.add_force(particles.forces.CohesiveSpamForce(p,p.nlists[1]))
    p.add_force(particles.forces.CollisionForce(p,p.nlists[0]))

    for nl in p.nlists:
        nl.build_nl_brute_force()
    cnt = 0


def update(dt):
    global cnt,fps,rebuild_nl 
    cnt += 1
    p.update()
    if cnt >= max_steps:
        pyglet.app.exit()
    pass

@s.win.event
def on_draw():
    s.fps = pyglet.clock.get_fps()
    s.redraw(p)    

@s.win.event
def on_key_press(symbol,modifiers):
    if symbol == pyglet.window.key.R:
        initialise()

def main():
    initialise()
    pyglet.clock.schedule_interval(update,1/20.0)
    pyglet.app.run()
    #loop()

main()
#profile.run('main()')


#Create some particles
#Use numerical intergration to step them forward in time
#The numerical integrator uses a subroutine that computes the forces
#SPHFORCES (subroutine)
## outputs 
#vdot 
#udot

## inputs
#Q
#P


