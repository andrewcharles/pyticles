#! /usr/local/bin/python

""" This is the main file for the python sph code
    We have a major memory leak!
    Write for 3 dimensions as much as possible
"""

import sys
#sys.path.append("/home/ac/pyglet-1.0")
import time
import particles
from pyglet import font
from pyglet import window
from pyglet import clock
from pyglet import app
import spyview
import profile

# Global variables
max_steps = 1000

p = particles.Particles()
nl_1 = particles.neighbour_list.NeighbourList(p,3.0)
nl_2 = particles.neighbour_list.NeighbourList(p,10.0)

#p.add_force(particles.forces.CollisionForce(p,p.nlists[0]))
#p.add_force(particles.forces.HookesForce(p,p.nlists[1]))

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
    nl_1 = particles.neighbour_list.NeighbourList(p,3.0)
    nl_2 = particles.neighbour_list.NeighbourList(p,10.0)
    p.nlists.append(nl_1)
    p.nlists.append(nl_2)
    p.nl_default = nl_2

#    p.add_force(particles.forces.SpamForce(p,p.nlists[1]))

    for nl in p.nlists:
        nl.build_nl_brute_force()
    cnt = 0


def update(dt):
    global cnt,fps,rebuild_nl 
    s.fps = clock.get_fps()
    cnt += 1
    p.update()
    if cnt >= max_steps:
        app.exit()
    pass

@s.win.event
def on_draw():
    s.redraw(p)    

@s.win.event
def on_key_press(symbol,modifiers):
    if symbol == window.key.R:
        initialise()

def main():
    initialise()
    clock.schedule_interval(update,1/60.0)
    app.run()
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


