#! /usr/local/bin/python

""" This is the main file for the particle system code.
    This does not use the SPH particle class currently - see spana
    for details of that. This one is for driving developments
    of the GUI and of fast generic particle system components.

    Copyright Andrew Charles 2008
    All rights reserved.
"""

import sys
import time
import particles
import pyglet
import _forces
from pyglet.window import mouse
import pview
import profile
from pyglet.gl import *
import controller
import gui

# Global variables
max_steps = 1000
NP1 = 2 
NP2 = 0 

p = particles.ParticleSystem(NP1,d=3,controllers=[controller.BodyForce()])
buttons = []
#b = particles.Particles(NP2)
#b.colour = 0.0,0.0,1.0

# nl_1 gives internal p interactions
# nl_2 gives p-b interactions
# there are no internal b interactions
#nl_1 = particles.neighbour_list.NeighbourList(p,3.0)
#nl_2 = particles.neighbour_list.NeighbourList(p,3.0,particle2=b)

s = pview.ParticleView()

cnt = 0
fps = 0
tstart = time.time()
rebuild_nl = 1

cmd_label = pyglet.text.Label("Command",font_name="Arial", \
            font_size=12,color =(240,0,220,244),x=500,y=460 )

act_label = pyglet.text.Label("Command",font_name="Arial", \
            font_size=12,color =(100,100,220,244),x=200,y=400 )

def update(dt):
    global cnt,fps,rebuild_nl 
    cnt += 1
    p.update()
    #b.update()
    if cnt >= max_steps:
        pyglet.app.exit()
    pass

@s.win.event
def on_draw():
    s.fps = pyglet.clock.get_fps()
    s.clear()
    s.redraw(p)
    cmd_label.draw()
    act_label.draw()
    for b in buttons:
        b.draw()    

@s.win.event
def on_key_press(symbol,modifiers):
    if symbol == pyglet.window.key.R:
        initialise()

@s.win.event
def on_mouse_motion(x,y,dx,dy):
    for b in buttons:
        if b.hit(x,y):
            #print "mouseover ",b.label
            cmd_label.text = b.label
            #cmd_label.x=b.x
            #cmd_label.y=b.y
    
@s.win.event
def on_mouse_press(x,y,button,modifiers):
    if button == mouse.LEFT:
        for b in buttons:
            if b.hit(x,y):
                b.activate()
                return
        p.create_particle(x/s.xmap,y/s.ymap)
        print "Creating particle at ",x/s.xmap,y/s.ymap


def clear_forces():
    global p 
    for nl in p.nlists:
        nl.forces=[]
        nl.nforce=0

def add_hookes():
    print "Adding spring force"
    p.nl_default.add_force(particles.forces.HookesForce(p,p,p.nl_default))

def add_grav():
    print "Adding gravity force"
    p.nl_default.add_force(particles.forces.Gravity(p,p,p.nl_default))

def create_ui():
    global buttons
    pyglet.resource.path.append('res')
    pyglet.resource.reindex()
    hookes_button = gui.Button()
    hookes_button.colour = 0.5,0.1,0.1
    hookes_button.x = 600
    hookes_button.y = 355
    hookes_button.label = "Springs"
    hookes_button.activate = add_hookes
    hookes_button.img = pyglet.resource.image('spring.png')

    grav_button = gui.Button()
    grav_button.colour = 0.1,0.5,0.1
    grav_button.x = 600 
    grav_button.y = 395
    grav_button.label = "grav"
    grav_button.activate = add_grav

    clear_button = gui.Button()
    clear_button.colour = 0.1,0.1,05.
    clear_button.x = 600
    clear_button.y = 435
    clear_button.activate=clear_forces
    clear_button.label = "clear"

    buttons.append(hookes_button)
    buttons.append(grav_button)
    buttons.append(clear_button)

def initialise():
    global p,nl_1,nl_2,cnt,buttons
    print "Restarting"
    p = particles.ParticleSystem(NP1,d=3
        ,controllers=[controller.BodyForce()]
        )
    
    nl_1 = particles.neighbour_list.NeighbourList(p,10.0)
    #nl_2 = particles.neighbour_list.NeighbourList(p,5.0)
    #nl_3 = particles.neighbour_list.NeighbourList(p,3.0,particle2=b)
    
    p.nlists.append(nl_1)
    #p.nlists.append(nl_2)
    #p.nlists.append(nl_3)
    p.nl_default = nl_1

    #nl_1.add_force(particles.forces.SpamForce(p,p.nlists[0]))
    #nl_1.add_force(particles.forces.CohesiveSpamForce(p,p.nlists[1]))
    #nl_1.add_force(particles.forces.CollisionForce(p,p,nl_1))
    #nl_1.add_force(particles.forces.HookesForce(p,p,nl_1))
    #nl_1.add_force(particles.forces.Gravity(p,p,nl_1))

    for nl in p.nlists:
        nl.build_nl_verlet()
    cnt = 0

    create_ui()

def main():
    #glEnable(GL_BLEND)
    initialise()
    pyglet.clock.schedule_interval(update,1/5.0)
    pyglet.app.run()
    #loop()

if __name__ == "__main__":
    main()


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


