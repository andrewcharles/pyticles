
# viewer module for particle system

import sys
sys.path.append("/home/ac/pyglet-1.0")
import particles
from pyglet import window
from pyglet import options
options['debug_gl'] = False
from pyglet.gl import *

from pylab import *
import pdb
from pyglet import font
from pyglet import image

sys.path.append("../../vasp")
import spkernel

WINDOW_WIDTH = 640
WINDOW_HEIGHT = 480

class ParticleView:
    " A particle viewer"
    def __init__(self,p):
        #config = Config(alpha_size=8)
        self.win = window.Window(WINDOW_WIDTH,WINDOW_HEIGHT,visible=False)
        self.ft = font.load('Arial',36)
        self.txt = font.Text(self.ft,"Hello pyglet")
        self.txt.color = (1.0, 0.0, 0.0, 1.0)
        self.win.set_visible()
        # multipliers to map system to screen coordinates
        self.xmap = WINDOW_WIDTH / particles.XMAX
        self.ymap = WINDOW_HEIGHT / particles.YMAX
        #self.img = image.load('p.png')  

    def draw_particles(self,p):
        """ issues the opengl commands to draw the 
            particle system """
        radius=5
        for i in range(p.n):
        #    self.img.blit(p.r[i,0],p.r[i,1])
            r = p.r[i,0] * self.xmap, p.r[i,1] * self.ymap
            glBegin(GL_POLYGON)
            glColor3f(1.0,0.0,0.0)
            for angle in range(6):
                a = radians(angle*60)
                glVertex2f(r[0]+sin(a)*radius,r[1]+cos(a)*radius)
            glEnd()
        
    def redraw(self,p,message):
        self.win.dispatch_events()
        glClear(GL_COLOR_BUFFER_BIT)
        
        self.win.clear()
        self.txt = font.Text(self.ft,message)
        self.txt.draw()
        
        #glLoadIdentity()
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
        self.draw_particles(p)
        self.win.flip()
