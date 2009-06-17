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
import forces
import pyglet
from pyglet.window import mouse
import pview
import profile
import neighbour_list
from pyglet.gl import *
import gui

# Global variables
MAX_STEPS = 10000
NP1 = 1
MAXN = 50

p = particles.SmoothParticleSystem(NP1,maxn=MAXN)
s = pview.ParticleView()

buttons = []
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
    p = particles.SmoothParticleSystem(NP1,maxn=MAXN)
    
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
    create_ui()


def add_hookes():
    print "Adding spring force"
    p.forces.append(forces.HookesForce(p,p,p.nl_default))

def add_grav():
    print "Adding gravity force"
    p.forces.append(forces.Gravity(p,p,p.nl_default))

def add_spam():
    print "Adding smooth particle hydrodynamic force"
    p.forces.append(forces.SpamForce(p,nl_1))

def add_spam_attract():
    print "Adding smooth particle hydrodynamic force"
    p.forces.append(forces.CohesiveSpamForce(p,nl_2))

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

    spam_button = gui.Button(
                    loc = (bx,300)
                    ,color = (0.9,0.3,0.5)
                    ,activate = add_spam
                    ,image = None
                    ,label = "button"
                )
    buttons.append(spam_button)


def clear_forces():
    global p 
    p.forces=[]

def update(dt):
    global cnt,fps,rebuild_nl 
    cnt += 1
    p.update()
    if cnt >= MAX_STEPS:
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



def main():
    initialise()
    pyglet.clock.schedule_interval(update,1/5.0)
    pyglet.app.run()

if __name__ == "__main__":
    main()




