"""
    Viewer module for particle system
    Copyright Andrew Charles 2008
    All rights reserved.

    To mess around with the various visualisation
    options edit the function redraw()

"""

import sys
import particles
import pyglet
pyglet.options['debug_gl'] = False
from pyglet.gl import *
import numpy

from pygarrayimage.arrayimage import ArrayInterfaceImage

from pylab import *
import pdb
import spkernel
import interpolate
import renderspam
import scipy
import scipy.ndimage
import pylab
from PIL import Image as pilmage
import scipy.misc.pilutil
from trackball_camera import TrackballCamera

WIDTH = 640
HEIGHT = 480
DEPTH = 480
RES = 2.1
width=WIDTH
height=HEIGHT


class ParticleView:
    " A particle viewer"
    def __init__(self):
        #config = Config(alpha_size=8)
        self.win = pyglet.window.Window(WIDTH,HEIGHT
                 ,visible=False,caption='Pyticles')

        self.sphere = gluNewQuadric()
        gluQuadricDrawStyle(self.sphere,GLU_LINE)
        glEnable(GL_DEPTH_TEST)
        glDisable(GL_CULL_FACE)

        # Set up labels. These labels are pretty tightly coupled to the
        # gui buttons
        self.fps = 0
        self.fpslab = pyglet.text.Label("fps label",font_name="Arial", \
            font_size=12,color =(255,0,0,255),x=10,y=4 )
        self.npart = pyglet.text.Label("np label",font_name="Arial", \
            font_size=12,color =(255,0,0,255),x=100,y=4 )
        self.nebs = pyglet.text.Label("ni label",font_name="Arial", \
            font_size=12,color =(255,0,0,255),x=10,y=18 )
       

        self.tbcam  = TrackballCamera(200.0) 
        self.win.on_resize = self.resize
        self.win.on_mouse_press = self.on_mouse_press
        self.win.on_mouse_drag = self.on_mouse_drag
        self.win.on_key_press = self.wsad
        

        self.win.set_visible()
        
        # multipliers to map system to screen coordinates
        self.xmap = WIDTH / 10.0 / float(particles.XMAX)
        self.ymap = HEIGHT / 10.0 / float(particles.YMAX)
        self.zmap = DEPTH / 10.0 / float(particles.ZMAX)
        

    def draw_box(self):
        x = float(particles.XMAX) * self.xmap
        y = float(particles.YMAX) * self.ymap
        z = float(particles.ZMAX) * self.zmap

        """
    
    (0,y,z) --------- (x,y,z)
       |\                | \
       | \               |  \
       |  \              |   \
       |   \             |    \
       |   (0,y,0) -------- (x,y,0)
       |    |            |      |
    (0,0,z) | ------- (x,0,z)   |
        \   |               \   |
         \  |                \  |
          \ |                 \ |
           \|                  \|
        (0,0,0) ----------- (x,0,0)

        """
        
        glBegin(GL_LINES)
        glColor3f(0.0,0.0,1.0)
        
        glVertex3f(0,0,0)
        glVertex3f(0,0,z)
        glVertex3f(0,y,z)
        glVertex3f(0,y,0)
        glVertex3f(0,0,0)
        glVertex3f(x,0,0)
        glVertex3f(x,y,0)
        glVertex3f(0,y,0)
        glVertex3f(0,y,z)
        glVertex3f(x,y,z)
        glVertex3f(x,0,z)
        glVertex3f(0,0,z)
        glVertex3f(0,0,0)

        glEnd()
        
        
    def draw_particles(self,p):
        """ issues the opengl commands to draw the 
            particle system """
        #glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
        radius=5
        for i in range(p.n):
            r = p.r[i,0] * self.xmap \
              , p.r[i,1] * self.ymap \
              , p.r[i,2] * self.zmap
            glColor3f(1.0, 1.0/(i+1), 0)
            glPushMatrix()
            glTranslatef(r[0],r[1],r[2])
            #glTranslatef(-1.,-1.,-1.)
            gluSphere(self.sphere,1.0,20,20)
            glPopMatrix()

            #glBegin(GL_POLYGON)
            #glColor3f(p.colour[0],p.colour[1],p.colour[2])
            #for angle in range(6):
            #    a = radians(angle*60)
            #    glVertex2f(r[0]+sin(a)*radius,r[1]+cos(a)*radius)
            #glEnd()

    def draw_neighbours(self,p):
        """ issues the opengl commands to draw 
            lines connecting neighbours """
        for i in range (p.nl_default.nip):
            glBegin(GL_LINES)
            glColor3f(0.3,0.3,0.0)
            a = p.nl_default.iap[i,0]
            b = p.nl_default.iap[i,1]
            glVertex2f(p.r[a,0]*self.xmap,p.r[a,1]*self.ymap)
            glVertex2f(p.r[b,0]*self.xmap,p.r[b,1]*self.ymap)
            glEnd()
  
    def hud(self,p):
        """ Draws the heads up display """
        k = 20
        xi = 10
        yi = 450

        self.set_ortho()
        glPushMatrix()
        glLoadIdentity()
        
        self.fpslab.text = "fps: %4.2f" %(self.fps)
        self.fpslab.draw()
        self.nebs.text = "n: "+str(p.n)
        self.nebs.draw()
        self.npart.text = "pairs: "+str(p.nl_default.nip)
        self.npart.draw()
        
        glPopMatrix()
        self.reset_perspective()

    def clear(self): 
        self.win.dispatch_events()
        glClear(GL_COLOR_BUFFER_BIT)
        glLoadIdentity()
        self.win.clear()
        
    def redraw(self,p):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        #self.draw_neighbours(p)
        self.draw_particles(p)
        self.draw_box()
        #self.hud(p)

    def set_ortho(self):
        """ Sets up an orthographic projection for drawing text and the like.
        """
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        #gluOrtho2D(0,width,0,height)
        glOrtho(0,width,0,height,-1,1)
        glMatrixMode(GL_MODELVIEW)

    def reset_perspective(self):
        """ Pops the projection matrix. Used for getting the
            3D perspective view back.
        """
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)

    def resize(self,width, height):
        """Setup 3D projection for window"""
        print self,width,height
        glViewport(0, 0, width,height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45, 1.0 * width/height, 0.001, 10000)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        self.tbcam.update_modelview()

    def norm1(self,x,maxx):
        """given x within [0,maxx], scale to a range [-1,1]."""
        return (2.0 * x - float(maxx)) / float(maxx)

    def on_mouse_press(self, x, y, button, modifiers):
        if button == pyglet.window.mouse.LEFT:
            self.tbcam.mouse_roll(
                self.norm1(x,width),
                self.norm1(y,height),
                False)
        elif button == pyglet.window.mouse.RIGHT:
            self.tbcam.mouse_zoom(
                self.norm1(x,width),
                self.norm1(y,height),
                False)

    def on_mouse_drag(self,x, y, dx, dy, buttons, modifiers):
        if buttons & pyglet.window.mouse.LEFT:
            if modifiers & pyglet.window.key.MOD_SHIFT:
                self.tbcam.mouse_roll(
                    self.norm1(x,width),
                    self.norm1(y,height))
            else: 
                self.tbcam.mouse_zoom(
                    self.norm1(x,width),
                    self.norm1(y,height))

    def wsad(self,symbol,modifiers):
        """ zoom and pan directly with the fps keys.
        """
        if symbol == pyglet.window.key.W:
            # zoom in
            print 'you called zoom_in'
        if symbol == pyglet.window.key.S:
            # zoom out
            print 'you called zoom_out'
        if symbol == pyglet.window.key.A:
            # pan left
            print 'you called pan left'
        if symbol == pyglet.window.key.D:
            # pan right
            print 'you called pan right'

