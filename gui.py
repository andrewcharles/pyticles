""" Graphical user interface elements using simple pyglet
    primitives.
"""

import pyglet
from pyglet.window import mouse
from pyglet.gl import *

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
