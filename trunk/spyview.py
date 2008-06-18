
# viewer module for particle system

import sys
import particles
import pyglet
pyglet.options['debug_gl'] = False
from pyglet.gl import *
import numpy

from pygarrayimage.arrayimage import ArrayInterfaceImage

from pylab import *
import pdb

sys.path.append("../../vasp")
import spkernel
import interpolate
import renderspam
import scipy
import scipy.ndimage
import pylab
#from PIL import Image as pilmage
# \todo install python image library
WINDOW_WIDTH = 640 #640
WINDOW_HEIGHT = 480 #480
RES =1
# res greater than 1 is still not registering in the right spot! 

class ParticleView:
    " A particle viewer"
    def __init__(self):
        #config = Config(alpha_size=8)
        self.win = pyglet.window.Window(WINDOW_WIDTH,WINDOW_HEIGHT,visible=False)
        self.fps = 0
        self.fpslab = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=10,y=400 )
        self.npart = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=10,y=380 )
        self.nebs = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=10,y=360 )
        self.maxvol = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=10,y=340 )
        
        self.win.set_visible()
        # multipliers to map system to screen coordinates
        self.xmap = WINDOW_WIDTH / float(particles.XMAX)
        self.ymap = WINDOW_HEIGHT / float(particles.YMAX)
        #self.img = image.load('p.png')
        self.fps = 0
        self.gridx,self.gridy = renderspam.get_grid_map(0.0,particles.XMAX,0.0,particles.YMAX,RES)  

    def render_density(self,p):
        """ create an image, and color the pixels,
            using the vasp module            
        """ 
        cutoff=5.0
        bounds = 0,particles.XMAX,0,particles.YMAX
        Z = interpolate.splash_map(self.gridx,self.gridy,p.r[0:p.n,:],p.m[0:p.n],p.rho[0:p.n],p.rho[0:p.n],p.h[0:p.n],bounds,cutoff)
        #Z=Z.transpose()
        #print rx
        #print ry
        #print Z.shape 
        # need to normalise wrt to whatever the color range is
        # and then figure out where the nan's are coming from
        #x,y = scipy.ogrid[0:particles.XMAX:RES,0:particles.YMAX:RES]
        rx,ry = scipy.mgrid[0:WINDOW_WIDTH:1,0:WINDOW_HEIGHT:1]
        #print x
        #print y
        #dx = x[1,0] - x[0,0]
        #dy = y[0,1] - y[0,0]
        #xvals = (rx - x[0,0])/dx
        #yvals = (ry - y[0,0])/dy
        #dx = self.gridx[1]-self.gridx[0]
        #dy = self.gridy[1]-self.gridy[0]
        #dx = WINDOW_WIDTH/self.gridx.shape[0]
        #dy = WINDOW_HEIGHT/self.gridy.shape[0]
        dx = float(RES*WINDOW_WIDTH/float(particles.XMAX))
        dy =float( RES*WINDOW_HEIGHT/float(particles.YMAX))
        xvals = (rx)/dx #(rx-RES?)
        yvals = (ry)/dy
        #yvals = yvals.transpose()
        G = numpy.array([xvals,yvals])
        #print "G"
        # map_coordinates seems to be creating some pixels
        # with values less than zero, which is fudging up
        # the rendering, if order >2 is used
        D = scipy.ndimage.map_coordinates(Z,G,order=1) 
        #print D.shape
        # need to create an image the size of the screen
        # and then interpolate based on our splashmap
        D = 255.*(D/amax(D))
        #D = D.transpose()
        D = numpy.flipud(D)
        #pylab.imshow(D)
        #pylab.colorbar()
        #pylab.show()
        D = numpy.cast['uint8'](D)
        #D = D.astype('uint8')
        D = scipy.ndimage.rotate(D,270,axes=(0,1))
        A = ArrayInterfaceImage(D)
        #img = A.texture
        img = A 
        img.blit(0,0,0)


    def draw_particles(self,p):
        """ issues the opengl commands to draw the 
            particle system """
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
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
        
        for i in range (p.nl_default.nip):
            """ draw lines connecting neighbours """
            glBegin(GL_LINES)
            glColor3f(1.0,1.0,1.0)
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
        self.fpslab.text = "fps: "+str(self.fps)
        self.fpslab.draw()
        self.nebs.text = "n: "+str(p.n)
        self.nebs.draw()
        self.npart.text = "pairs: "+str(p.nl_default.nip)
        self.npart.draw()
        self.maxvol.text = "V_max: "+ str(max(p.m[0:p.n]/p.rho[0:p.n]))
        self.maxvol.draw()
 
        
    def redraw(self,p):
        self.win.dispatch_events()
        glClear(GL_COLOR_BUFFER_BIT)
        self.win.clear()
        #self.render_density(p)
        self.hud(p)
        self.draw_particles(p)
