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
NP1 = 20
NP2 = 3

p = particles.Particles(NP1)
b = particles.Particles(NP2)
b.colour = 0.0,0.0,1.0
# nl_1 gives internal p interactions
# nl_2 gives p-b interactions
# there are no internal b interactions
nl_1 = particles.neighbour_list.NeighbourList(p,3.0)
nl_2 = particles.neighbour_list.NeighbourList(p,3.0,particle2=b)

s = spyview.ParticleView()

cnt = 0
fps = 0
tstart = time.time()
rebuild_nl = 1

def initialise():
    global p,nl_1,nl_2,cnt
    print "Restarting"
    p = particles.Particles(NP1)
    b = particles.Particles(NP2)
    
    nl_1 = particles.neighbour_list.NeighbourList(p,3.0)
    nl_2 = particles.neighbour_list.NeighbourList(p,3.0,particle2=b)
    
    p.nlists.append(nl_1)
    p.nlists.append(nl_2)
    p.nl_default = nl_1

    #p.add_force(particles.forces.SpamForce(p,p.nlists[0]))
    #p.add_force(particles.forces.CohesiveSpamForce(p,p.nlists[1]))
    nl_1.add_force(particles.forces.CollisionForce(p,p,p.nlists[0]))
    nl_2.add_force(particles.forces.HookesForce(p,b,p.nlists[1]))

    for nl in p.nlists:
        nl.build_nl_verlet()
    cnt = 0


def update(dt):
    global cnt,fps,rebuild_nl 
    cnt += 1
    p.update()
    b.update()
    if cnt >= max_steps:
        pyglet.app.exit()
    pass

@s.win.event
def on_draw():
    s.fps = pyglet.clock.get_fps()
    s.clear()
    s.redraw(b)    
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


