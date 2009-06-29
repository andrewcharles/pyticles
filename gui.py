""" Graphical user interface elements using simple pyglet
    primitives.
"""

import pyglet
from pyglet.window import mouse
from pyglet.gl import *

def nothing():
    pass

class Button:
    """ A button.

        __init__() -- create a button with default values

    """
    
    def __init__(self
                ,loc = (30,30)
                ,color = (0.5,0.5,0.5)
                ,mcolor = (0.9,0.9,0.9)
                ,activate = nothing()
                ,image = None
                ,label = "button"
                ):
        self.color = color
        self.mcolor = mcolor
        self.x = loc[0]
        self.y = loc[1]
        self.width = 32
        self.height = 32
        self.activate = activate
        self.img = image
        self.label=label

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
        if self.img:
            self.img.blit(self.x,self.y)
        else:
            glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
            glBegin(GL_POLYGON)
            glColor3f(self.color[0],self.color[1],self.color[2])
            glVertex2f(self.x,self.y)
            glVertex2f(self.x,self.y+self.height)
            glVertex2f(self.x+self.width,self.y+self.height)
            glVertex2f(self.x+self.width,self.y)
            glVertex2f(self.x,self.y)
            glEnd()


