#! /usr/local/bin/python
""" A test case to show we is leaking memory.
"""
import pyglet
from pyglet.gl import *
from pygarrayimage.arrayimage import ArrayInterfaceImage
import numpy
import scipy.ndimage

WINDOW_WIDTH =640 
WINDOW_HEIGHT =480 
win = pyglet.window.Window(WINDOW_WIDTH,WINDOW_HEIGHT,visible=False)

def redraw():
    glClear(GL_COLOR_BUFFER_BIT)
    D = numpy.random.random([640,480])
    D = D*255
    D = numpy.cast['uint8'](D)
    D = scipy.ndimage.rotate(D,270,axes=(0,1))
    A = ArrayInterfaceImage(D)
    img = A.image_data 
    img.blit(0,0)

def update(dt):
    pass #the salt

@win.event
def on_draw():
    redraw()    

def main():
    print pyglet.version
    win.set_visible()
    pyglet.clock.schedule_interval(update,1/60.0)
    pyglet.app.run()

main()

