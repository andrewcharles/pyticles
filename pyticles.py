#! /usr/local/bin/python

""" This is the main file for the python sph code
    Copyright Andrew Charles 2008
    All rights reserved.
"""

import sys
#sys.path.append("/home/ac/pyglet-1.0")
import time
import particles
import pyglet
from pyglet.window import mouse
import spyview
import profile
from pyglet.gl import *

# Global variables
max_steps = 1000
NP1 = 2 
NP2 = 0 

p = particles.Particles(NP1)
buttons = []
#b = particles.Particles(NP2)
#b.colour = 0.0,0.0,1.0

# nl_1 gives internal p interactions
# nl_2 gives p-b interactions
# there are no internal b interactions
#nl_1 = particles.neighbour_list.NeighbourList(p,3.0)
#nl_2 = particles.neighbour_list.NeighbourList(p,3.0,particle2=b)

s = spyview.ParticleView()

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

class Button:
    """ A button.
    """
    
    def __init__(self):
        self.colour=0.1,0.1,0.1
        self.mcolour=0.9,0.9,0.9
        self.x = 2
        self.y = 2
        self.width = 40
        self.height = 40
        self.img = "none"

    def activate(self):
            print "Hit the button"

    def hit(self,x,y):
        """ Checks if the given coords are in the
            button's hitbox.
        """
        if( (x>self.x) and (y>self.y)
        and (x<self.x+self.width) and (y<self.y+self.width)):
            #self.activate()
            return True
    
    def draw(self):
        if self.img != "none":
            self.img.blit(self.x,self.y)
        else:
            glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
            glBegin(GL_POLYGON)
            glColor3f(self.colour[0],self.colour[1],self.colour[2])
            glVertex2f(self.x,self.y)
            glVertex2f(self.x,self.y+self.height)
            glVertex2f(self.x+self.width,self.y+self.height)
            glVertex2f(self.x+self.width,self.y)
            glVertex2f(self.x,self.y)
            glEnd()


def clear_forces():
    global p 
    for nl in p.nlists:
        nl.forces=[]
        nl.nforce=0

def add_hookes():
    print "Adding spring force"
    p.nl_default.add_force(particles.forces.HookesForce(p,p,p.nl_default))

def add_spam():
    print "Adding spam force"
    p.nl_default.add_force(particles.forces.SpamForce(p,p.nl_default))

def create_ui():
    global buttons
    pyglet.resource.path.append('res')
    pyglet.resource.reindex()
    hookes_button = Button()
    hookes_button.colour = 0.5,0.1,0.1
    hookes_button.x = 600
    hookes_button.y = 355
    hookes_button.label = "Springs"
    hookes_button.activate = add_hookes
    hookes_button.img = pyglet.resource.image('spring.png')

    sph_button = Button()
    sph_button.colour = 0.1,0.5,0.1
    sph_button.x = 600 
    sph_button.y = 395
    sph_button.label = "Spam"
    sph_button.activate = add_spam

    clear_button = Button()
    clear_button.colour = 0.1,0.1,05.
    clear_button.x = 600
    clear_button.y = 435
    clear_button.activate=clear_forces
    clear_button.label = "clear"

    buttons.append(hookes_button)
    buttons.append(sph_button)
    buttons.append(clear_button)

def initialise():
    global p,nl_1,nl_2,cnt,buttons
    print "Restarting"
    p = particles.Particles(NP1)
    
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
    nl_1.add_force(particles.forces.Gravity(p,p,nl_1))

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


