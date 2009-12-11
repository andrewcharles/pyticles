#! /usr/local/bin/python

""" 
    A front end for the python SPH code.
    

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
from controller import BodyForce

# Global variables
MAX_STEPS = 10000
NP1 = 27
MAXN = 50
XMAX = 20 #64
YMAX = 20 #48
ZMAX = 20
VMAX = 0.0
dt = 0.1
SPACING = 3.0

# Model
p = particles.SmoothParticleSystem(NP1,maxn=MAXN,rinit='grid',side=(3,3,3)
    ,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX)

# Views
s = pview.SmoothParticleView(p)
#info = pview.ParticleInfo(p)

buttons = []
cnt = 0
fps = 0
tstart = time.time()
rebuild_nl = 1
bg_color = (0.3,0.3,0.3,1.0)

act_label = pyglet.text.Label("Command",font_name="Arial", \
            font_size=12,color =(100,100,220,244),x=200,y=20 )

default_group = pyglet.graphics.Group()

def initialise():
    global p,nl_1,nl_2,cnt,buttons
    print "Restarting"
    p = particles.SmoothParticleSystem(NP1,maxn=MAXN,rinit='grid',side=(3,3,3)
        ,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX,vmax=VMAX)
    nl_1 = neighbour_list.VerletList(p,cutoff=5.0)
    nl_2 = neighbour_list.VerletList(p,cutoff=10.0)
    p.nlists.append(nl_1)
    p.nlists.append(nl_2)
    p.nl_default = nl_1

    p.forces.append(forces.SpamForce(p,nl_1))
    p.forces.append(forces.CohesiveSpamForce(p,nl_2))
    #p.controllers.append(BodyForce(p,nl_2))

    for nl in p.nlists:
        nl.build()

    cnt = 0
    create_ui()

# Define GUI hook functions

def add_hookes():
    print "Adding spring force"
    p.forces.append(forces.HookesForce(p,p.nl_default))
    nl_1.build()

def add_grav():
    print "Adding gravity force"
    p.forces.append(forces.Gravity(p,nl_1))
    nl_1.build()

def add_spam():
    print "Adding smooth particle hydrodynamic force"
    p.forces.append(forces.SpamForce(p,nl_1))
    nl_1.build()

def add_spam_attract():
    print "Adding smooth particle hydrodynamic force"
    p.forces.append(forces.CohesiveSpamForce(p,nl_2))
    nl_2.build()

def add_particle():
    p.create_particle((0.0,YMAX,0.0),v=(0.1,0.0,0.0))
    nl_1.build()
    nl_2.build()

def inc_dt():
    global dt
    dt *= 1.2
    p.dt = dt

def dec_dt():
    global dt
    dt *= 0.8
    p.dt = dt

def create_ui():
    global buttons
    pyglet.resource.path.append('res')
    pyglet.resource.reindex()

    # GUI parameters are here for now
    # Assume a little 640 by 480 window
    # Put the buttons on the bottom like an
    # MMO action bar. This is nice and simple
    # so we don't need a layot manager.
    # 
    # (abar_x,abar_y) - action bar bottom left corner
    #  button_width
    # -------------------------------
    button_height = 32
    button_width = 64
    bx = [580,580,580,580,580,580,580]
    by = [300,260,220,180,140,100,60,20]

    # add button
    add_button = gui.Button(
                    loc =(bx[0],by[0]),
                    size = (button_width,button_height),
                    color = (0.5,0.5,0.5),
                    activate = add_particle,
                    labeltext = "Add",
                    group = default_group
                    #image = pyglet.resource.image('clear_button.png')
                    )
    buttons.append(add_button)

    # clear forces button
    clear_button = gui.Button(
                    loc =(bx[1],by[1]),
                    size = (button_width,button_height),
                    color = (0.5,0.5,0.5),
                    activate=clear_forces,
                    labeltext = "Clear",
                    group = default_group
                    #image = pyglet.resource.image('clear_button.png')
                    )
    buttons.append(clear_button)

    # spam button
    spam_button = gui.Button(
                    loc = (bx[2],by[2]),
                    size = (button_width,button_height),
                    color = (0.5,0.5,0.5),
                    activate = add_spam,
                    labeltext = 'Spam'
                    #image = pyglet.resource.image('spam_button_yellow.png')
                  )
    buttons.append(spam_button)

    # increase dt
    dt_up = gui.Button(
                loc = (bx[3],by[3]),
                size = (button_width,button_height),
                color = (0.5,0.5,0.5),
                activate = inc_dt,
                image = None,
                labeltext = "Faster"
                )
    buttons.append(dt_up)

    # decrease dt
    dt_down = gui.Button(
                loc = (bx[4],by[4]),
                size = (button_width,button_height),
                color = (0.5,0.5,0.5),
                activate = dec_dt,
                image = None,
                labeltext = "Slower"
                )
    buttons.append(dt_down)

    # decrease dt
    dt_down = gui.Button(
                loc = (bx[5],by[5]),
                size = (button_width,button_height),
                color = (0.5,0.5,0.5),
                activate = None,
                image = None,
                labeltext = "Blank"
                )
    buttons.append(dt_down)

    # decrease dt
    dt_down = gui.Button(
                loc = (bx[6],by[6]),
                size = (button_width,button_height),
                color = (0.5,0.5,0.5),
                activate = None,
                image = None,
                labeltext = "Blank"
                )
    buttons.append(dt_down)

def clear_forces():
    global p 
    p.forces=[]

def update(t):
    global cnt,fps,rebuild_nl,dt 
    cnt += 1
    p.update(dt)
    if cnt >= MAX_STEPS:
        pyglet.app.exit()
    p.steps += 1
    pass

@s.win.event
def on_draw():
    s.fps = pyglet.clock.get_fps()
    s.clear()
    s.redraw(p)
        
    # Draw UI
    s.set_ortho()
    glPushMatrix()
    glLoadIdentity()
    act_label.draw()
    for b in buttons:
        b.draw()    
    
    glFlush()
    glPopMatrix()
    s.reset_perspective()
    glClearColor(*bg_color)

@s.win.event
def on_key_press(symbol,modifiers):
    if symbol == pyglet.window.key.R:
        initialise()

@s.win.event
def on_mouse_motion(x,y,dx,dy):
    for b in buttons:
        if b.hit(x,y):
            act_label.text = b.label.text
    
@s.win.event
def on_mouse_press(x,y,button,modifiers):
    """ If we hit a button, activate it. Otherwise, we have
        hit somewhere in space, so create a particle.
        Will need the concept of depth when the display is 3D.
    """
    if button == mouse.LEFT:
        for b in buttons:
            if b.hit(x,y):
                b.activate()
                return
        #p.create_particle(x/s.xmap,y/s.ymap)
        #for nl in p.nlists:
         #   nl.build()
        #print "Creating particle at ",x/s.xmap,y/s.ymap

def main():
    initialise()
    pyglet.clock.schedule_interval(update,1/10.0)
    pyglet.clock.schedule_interval(s.update_eye,1/2.0)
    pyglet.app.run()

if __name__ == "__main__":
    main()




