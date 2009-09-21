#! /usr/local/bin/python

""" This is the main file for the particle system code.
    This does not use the SPH particle class currently - see spana
    for details of that. This one is for driving developments
    of the GUI and of fast generic particle system components.

    Copyright Andrew Charles 2009
    All rights reserved.

    Global variables
    ----------------
    max_steps -- the program will end when this number of steps
        have been taken.
    timestep_size -- length of simulated timestep

"""

import sys
import time
import particles
import pyglet
import forces
#import c_forces as forces
sys.path.append('/Users/acharles/masters/active/fsph')
#import hello
#import collision
from pyglet.window import mouse
import pview
import profile
from pyglet.gl import *
import controller
import gui

# Global variables
max_steps = 100000
timestep_size = 0.01
NP1 = 9 
NP2 = 0 

p = particles.ParticleSystem(NP1,d=3,controllers=[controller.BodyForce()])
buttons = []

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
    p.forces.append(forces.Gravity(p,nl_1))

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
    #s.clear()
    s.redraw(p)
    draw_ui()

@s.win.event
def on_key_press(symbol,modifiers):
    if symbol == pyglet.window.key.R:
        initialise()

@s.win.event
def on_mouse_motion(x,y,dx,dy):
    for b in buttons:
        if b.hit(x,y):
            #print "mouseover ",b.label
            cmd_label.text = b.label.text
            cmd_label.x=b.x
            cmd_label.y=b.y
    
#@s.win.event
#def on_mouse_press(x,y,button,modifiers):
#    if button == mouse.LEFT:
#        for b in buttons:
#            if b.hit(x,y):
#                b.activate()
#                return
#        p.create_particle(x/s.xmap,y/s.ymap)
#        print "Creating particle at ",x/s.xmap,y/s.ymap
#        p.nlists[0].build()

def clear_forces():
    global p 
    p.forces = []
    p.controllers = []

def add_hookes():
    print "Adding spring force."
    p.forces.append(forces.HookesForce(p,p.nl_default))
    p.nlists[0].build()

def add_grav():
    print "Adding gravity force."
    p.forces.append(forces.Gravity(p,p.nl_default))
    p.nlists[0].build()

def add_collisions():
    print "Adding collision force."
    p.forces.append(forces.CollisionForce(p,p.nl_default))
    p.nlists[0].build()

def add_body_force():
    c = controller.BodyForce() 
    p.controllers.append(c)
    c.bind_particles(p)
    p.nlists[0].build()

def draw_ui():
    s.set_ortho()
    glPushMatrix()
    glLoadIdentity()
    cmd_label.draw()
    act_label.draw()
    for b in buttons:
        b.draw()    
    glFlush()
    glPopMatrix()
    s.reset_perspective()

def create_ui():
    global buttons
    pyglet.resource.path.append('res')
    pyglet.resource.reindex()
    bx = 600
    # spring button
    hookes_button = gui.Button(
        labeltext="Springs",
        color = (0.5,0.1,0.1),
        loc = (bx,355),
        activate = add_hookes,
        image = pyglet.resource.image('spring.png')
        )
    buttons.append(hookes_button)

    # gravity button
    grav_button = gui.Button(
        color = (0.1,0.5,0.1),
        loc = (bx,395),
        labeltext = "grav",
        activate = add_grav
        )
    buttons.append(grav_button)

    # clear forces button
    clear_button = gui.Button(
        color = (0.1,0.1,0.5),
        loc = (bx,435),
        activate = clear_forces,
        labeltext = "clear"
        )
    buttons.append(clear_button)

    collision_button = gui.Button(
                             loc = (bx,300)
                            ,color = (0.9,0.3,0.5)
                            ,activate = add_collisions
                            ,image = pyglet.resource.image('collision.png')
                            ,labeltext = "collision")
    buttons.append(collision_button)

    body_button = gui.Button(
                             loc = (bx,200)
                            ,color = (0.1,0.3,0.5)
                            ,activate = add_body_force
                            ,image = None
                            ,labeltext = "body force")
    buttons.append(body_button)

    # increase dt
    dt_up = gui.Button(
                loc = (bx,150)
                ,color = (1.0,0.0,0.0)
                ,activate = inc_dt
                ,image = None
                ,labeltext = "dt_up"
                )
    buttons.append(dt_up)

    # decrease dt
    dt_down = gui.Button(
                loc = (bx,100)
                ,color = (0.0,0.0,1.0)
                ,activate = dec_dt
                ,image = None
                ,labeltext = "dt_down"
                )
    buttons.append(dt_down)

def inc_dt():
    global dt
    dt *= 1.2

def dec_dt():
    global dt
    dt *= 0.8

def main():
    # We use the pyglet event loop to drive the application
    initialise()
    #pyglet.clock.schedule(update)
    pyglet.clock.schedule_interval(update,1/5.0)
    pyglet.app.run()

if __name__ == "__main__":
    main()


# CODE FOR MULTIPLE INTERACTING PARTICLE SYSTEMS
#b = particles.Particles(NP2)
#b.colour = 0.0,0.0,1.0
# nl_1 gives internal p interactions
# nl_2 gives p-b interactions
# there are no internal b interactions
#nl_1 = particles.neighbour_list.NeighbourList(p,3.0)
#nl_2 = particles.neighbour_list.NeighbourList(p,3.0,particle2=b)
