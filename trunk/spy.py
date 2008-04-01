#! /usr/local/bin/python

# This is the main file for the python sph code
# Write for 3 dimensions as much as possible

# We have a major memory leak!

import sys
#sys.path.append("/home/ac/pyglet-1.0")
import time
import particles
from pyglet import font
from pyglet import window
from pyglet import clock
import spyview

# Global variables

p = particles.Particles()

# a long and a short neighbour list
nl_1 = particles.neighbour_list.NeighbourList(p,3.0)
nl_2 = particles.neighbour_list.NeighbourList(p,10.0)

p.nlists.append(nl_1)
p.nlists.append(nl_2)
p.nl_default = nl_2

for nl in p.nlists:
    nl.build_nl_brute_force()

p.add_force(particles.forces.CollisionForce(p,p.nlists[0]))
p.add_force(particles.forces.HookesForce(p,p.nlists[1]))
p.add_force(particles.forces.SpamForce(p,p.nlists[1]))

s = spyview.ParticleView(p)

cnt = 0
fps = 0
tstart = time.time()
rebuild_nl = 1

def initialise():
    print "wtf"

# MAIN CONTROL LOOP
# should just use the Pyglet main loop
# step
def update():
    global cnt,fps,rebuild_nl 
    dt = clock.tick()
    fps = clock.get_fps()
    cnt += 1
    p.update()
    s.redraw(p,str(fps))    
    return True

def on_key_press(symbol,modifiers):
    if symbol == window.key.R:
        initialise()

def loop():
    while not s.win.has_exit:
    #for t in range(0,500):
        update()

def main():
    loop()

main()

# what integrator

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


