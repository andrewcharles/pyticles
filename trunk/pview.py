"""
    Viewer module for particle system
    Copyright Andrew Charles 2008
    All rights reserved.
    This module is new BSD licensed.

    To mess around with the various (all 3...) visualisation
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
WINDOW_WIDTH = 640
WINDOW_HEIGHT = 480
RES = 2.1 


class ParticleView:
    " A particle viewer"
    def __init__(self):
        #config = Config(alpha_size=8)
        self.win = pyglet.window.Window(WINDOW_WIDTH,WINDOW_HEIGHT
                 ,visible=False,caption='Pyticles')
        self.fps = 0
        self.fpslab = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=10,y=4 )
        self.npart = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=100,y=4 )
        self.nebs = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=10,y=18 )
        self.maxvol = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=100,y=18 )
        
        self.win.set_visible()
        # multipliers to map system to screen coordinates
        self.xmap = WINDOW_WIDTH / float(particles.XMAX)
        self.ymap = WINDOW_HEIGHT / float(particles.YMAX)
        self.fps = 0
        # gridx,gridy is a grid at the rendering resolution
        self.gridx,self.gridy = renderspam.get_grid_map(0.0,particles.XMAX, \
        0.0,particles.YMAX,RES)
        # rx,ry is a pixel grid
        rx,ry = scipy.mgrid[0:WINDOW_WIDTH:1,0:WINDOW_HEIGHT:1]
        dx = RES*WINDOW_WIDTH/float(particles.XMAX)
        dy = RES*WINDOW_HEIGHT/float(particles.YMAX)
        #xvals = (rx-dx/RES)/dx #(rx-RES?)
        xvals = (rx-dx/2)/dx #(rx-RES?)
        yvals = (ry-dy/2)/dy
        self.G = numpy.array([xvals,yvals])
          

    def draw_particles(self,p):
        """ issues the opengl commands to draw the 
            particle system """
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
        radius=5
        for i in range(p.n):
            r = p.r[i,0] * self.xmap, p.r[i,1] * self.ymap
            glBegin(GL_POLYGON)
            glColor3f(p.colour[0],p.colour[1],p.colour[2])
            for angle in range(6):
                a = radians(angle*60)
                glVertex2f(r[0]+sin(a)*radius,r[1]+cos(a)*radius)
            glEnd()

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
        self.fpslab.text = "fps: %4.2f" %(self.fps)
        self.fpslab.draw()
        self.nebs.text = "n: "+str(p.n)
        self.nebs.draw()
        self.npart.text = "pairs: "+str(p.nl_default.nip)
        self.npart.draw()

    def clear(self): 
        self.win.dispatch_events()
        glClear(GL_COLOR_BUFFER_BIT)
        glLoadIdentity()
        self.win.clear()
        
    def redraw(self,p):
        self.hud(p)
        self.draw_neighbours(p)
        self.draw_particles(p)



class SmoothParticleView(ParticleView):
    " A particle viewer"
    def __init__(self):
        #config = Config(alpha_size=8)
        self.win = pyglet.window.Window(WINDOW_WIDTH,WINDOW_HEIGHT,visible=False,caption='Pyticles')
        self.fps = 0
        self.fpslab = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=10,y=4 )
        self.npart = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=100,y=4 )
        self.nebs = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=10,y=18 )
        self.maxvol = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=100,y=18 )
        
        self.win.set_visible()
        # multipliers to map system to screen coordinates
        self.xmap = WINDOW_WIDTH / float(particles.XMAX)
        self.ymap = WINDOW_HEIGHT / float(particles.YMAX)
        #self.img = pyglet.image.load('p.png')
        self.fps = 0
        # gridx,gridy is a grid at the rendering resolution
        self.gridx,self.gridy = renderspam.get_grid_map(0.0,particles.XMAX, \
        0.0,particles.YMAX,RES)
        # rx,ry is a pixel grid
        rx,ry = scipy.mgrid[0:WINDOW_WIDTH:1,0:WINDOW_HEIGHT:1]
        dx = RES*WINDOW_WIDTH/float(particles.XMAX)
        dy = RES*WINDOW_HEIGHT/float(particles.YMAX)
        #xvals = (rx-dx/RES)/dx #(rx-RES?)
        xvals = (rx-dx/2)/dx #(rx-RES?)
        yvals = (ry-dy/2)/dy
        self.G = numpy.array([xvals,yvals])
          

    def render_density(self,p):
        """ Generates a splashmap, then creates an image the
            size of the viewing window and interpolates the splashmap
            onto it.
            Currently dependent on vasp modules. 
        """ 
        glColor3f(1.0,1.0,1.0)
        #cutoff=p.nl_default.cutoff_radius
        cutoff=5
        bounds = 0,particles.XMAX,0,particles.YMAX
        Z = interpolate.splash_map(self.gridx,self.gridy,p.r[0:p.n,:], \
            p.m[0:p.n],p.rho[0:p.n],p.rho[0:p.n],p.h[0:p.n],bounds,cutoff)
        D = scipy.ndimage.map_coordinates(Z,self.G,order=1,mode='nearest', \
            prefilter=False)
        #D = numpy.random.random([640,480])
        # need to create an image the size of the screen
        # and then interpolate based on our splashmap
        D = numpy.rot90(D)
        #print amax(D/amax(D))
        D = 255.*(D/amax(D))
        D = numpy.flipud(D)
        D = numpy.cast['uint8'](D)
        A = ArrayInterfaceImage(D)
        img = A.image_data
        #A = pilmage.fromarray(D,'L')
        #img = pyglet.image.load('wtf.png',A)
        img.blit(0,0)

    def draw_particles(self,p):
        """ issues the opengl commands to draw the 
            particle system """
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
        radius=5
        for i in range(p.n):
        #    self.img.blit(p.r[i,0],p.r[i,1])
            r = p.r[i,0] * self.xmap, p.r[i,1] * self.ymap
            glBegin(GL_POLYGON)
            glColor3f(p.colour[0],p.colour[1],p.colour[2])
            for angle in range(6):
                a = radians(angle*60)
                glVertex2f(r[0]+sin(a)*radius,r[1]+cos(a)*radius)
            glEnd()

    def draw_particle_sprites(self,p):
        """ Draws particles as sprites.
        """
        for i in range(p.n):
            self.img.blit(self.xmap*p.r[i,0],self.ymap*p.r[i,1])

        
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
        self.fpslab.text = "fps: %4.2f" %(self.fps)
        self.fpslab.draw()
        self.nebs.text = "n: "+str(p.n)
        self.nebs.draw()
        self.npart.text = "pairs: "+str(p.nl_default.nip)
        self.npart.draw()
        self.maxvol.text = "max vol: %4.2f" %(max(p.m[0:p.n]/p.rho[0:p.n]))
        self.maxvol.draw()

    def clear(self): 
        self.win.dispatch_events()
        glClear(GL_COLOR_BUFFER_BIT)
        glLoadIdentity()
        self.win.clear()
        
    def redraw(self,p):
        #self.win.dispatch_events()
        #glClear(GL_COLOR_BUFFER_BIT)
        #glLoadIdentity()
        #self.render_density(p)
        self.hud(p)
        #self.draw_neighbours(p)
        self.draw_particles(p)
        #self.draw_particle_sprites(p)

