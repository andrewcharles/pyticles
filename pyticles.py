#! /usr/local/bin/python

""" This is the main file for the particle system code.
    This does not use the SPH particle class currently - see spana
    for details of that. This one is for driving developments
    of the GUI and of fast generic particle system components.

    Create some particles
    Use numerical intergration to step them forward in time
    The numerical integrator uses a subroutine that computes the forces
    SPHFORCES (subroutine)
    ## outputs 
    vdot 
    udot
    Q
    P

    Copyright Andrew Charles 2008
    All rights reserved.
"""

import sys
import time
import particles
import pyglet
import forces
#import c_forces as forces
from pyglet.window import mouse
import pview
import profile
from pyglet.gl import *
import controller
import gui

# Global variables
max_steps = 1000
timestep_size = 0.01
NP1 = 2 
NP2 = 0 

p = particles.ParticleSystem(NP1,d=3,controllers=[controller.BodyForce()])
buttons = []

# CODE FOR MULTIPLE INTERACTING PARTICLE SYSTEMS
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


def initialise():
    global p,nl_1,nl_2,cnt,buttons
    print "Restarting"
    p = particles.ParticleSystem(NP1,d=3
        ,controllers=[controller.BodyForce()]
        )
    
    nl_1 = particles.neighbour_list.NeighbourList(p)
    p.nlists.append(nl_1)
    p.nl_default = nl_1
    
    p.forces.append(forces.CollisionForce(p,nl_1))
    #p.forces.append(forces.HookesForce(p,nl_1))
    #p.forces.append(forces.Gravity(p,nl_1))

    for nl in p.nlists:
        nl.build()
    cnt = 0

    create_ui()


def update(dt):
    global cnt,fps,rebuild_nl 
    cnt += 1
    p.update(timestep_size)
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
    p.forces = []
    p.controllers = []

def add_hookes():
    print "Adding spring force"
    p.forces.append(forces.HookesForce(p,p.nlists[0]))
    p.nlists[0].build()

def add_grav():
    print "Adding gravity force"
    p.forces.append(forces.Gravity(p,p.nl_default))
    p.nlists[0].build()

def add_collisions():
    p.forces.append(forces.CollisionForce(p,p.nl_default))
    p.nlists[0].build()

def add_body_force():
    c = controller.BodyForce() 
    p.controllers.append(c)
    c.bind_particles(p)
    p.nlists[0].build()


def create_ui():
    global buttons
    pyglet.resource.path.append('res')
    pyglet.resource.reindex()

    bx = 600

    # spring button
    hookes_button = gui.Button()
    hookes_button.colour = 0.5,0.1,0.1
    hookes_button.x = bx
    hookes_button.y = 355
    hookes_button.label = "Springs"
    hookes_button.activate = add_hookes
    hookes_button.img = pyglet.resource.image('spring.png')
    buttons.append(hookes_button)

    # gravity button
    grav_button = gui.Button()
    grav_button.colour = 0.1,0.5,0.1
    grav_button.x = bx
    grav_button.y = 395
    grav_button.label = "grav"
    grav_button.activate = add_grav
    buttons.append(grav_button)

    # clear forces button
    clear_button = gui.Button()
    clear_button.color = 0.1,0.1,05.
    clear_button.x = bx
    clear_button.y = 435
    clear_button.activate=clear_forces
    clear_button.label = "clear"
    buttons.append(clear_button)

    collision_button = gui.Button(
                             loc = (bx,300)
                            ,color = (0.9,0.3,0.5)
                            ,activate = add_collisions
                            ,image = "collision.png"
                            ,label = "collision")
    buttons.append(collision_button)

    body_button = gui.Button(
                             loc = (bx,200)
                            ,color = (0.1,0.3,0.5)
                            ,activate = add_body_force
                            ,image = None
                            ,label = "button")
    buttons.append(body_button)

def main():
    #glEnable(GL_BLEND)
    initialise()
    pyglet.clock.schedule_interval(update,1/5.0)
    pyglet.app.run()
    #loop()

if __name__ == "__main__":
    main()


