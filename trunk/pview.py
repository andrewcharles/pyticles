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

WINDOW_WIDTH = 640
WINDOW_HEIGHT = 480
WINDOW_DEPTH = 480
RES = 2.1
width=WINDOW_WIDTH
height=WINDOW_HEIGHT


class ParticleView:
    " A particle viewer"
    def __init__(self):
        #config = Config(alpha_size=8)
        self.win = pyglet.window.Window(WINDOW_WIDTH,WINDOW_HEIGHT
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
        self.win.on_key_press = self.wsad_press
        
        self.win.set_visible()
        
        # multipliers to map system to screen coordinates
        self.xmap = WINDOW_WIDTH / float(particles.XMAX)
        self.ymap = WINDOW_HEIGHT / float(particles.YMAX)
        self.zmap = WINDOW_DEPTH / float(particles.ZMAX)

        # Where is the center of the rectangle?
        cx = float(particles.XMAX/2.0) * self.xmap
        cy = float(particles.YMAX/2.0) * self.ymap
        cz = float(particles.ZMAX/2.0) * self.zmap
        self.tbcam.cam_focus = (cx,cy,cz)

        # zoom out
        self.tbcam.cam_eye[2] = -600
        self.tbcam.cam_eye[1] = -0
        self.tbcam.cam_eye[0] = -0
        self.tbcam.update_modelview()
        

    def draw_box(self):
        x = WINDOW_WIDTH
        y = WINDOW_HEIGHT
        z = WINDOW_DEPTH

        """
    
    (0,y,z) --------- (x,y,z)
       |\                | \
       | \          3    |  \
       |  \              |   \
       |   \             |    \
       |   (0,y,0) -------- (x,y,0)
       | 1  |            |      |
    (0,0,z) | ------- (x,0,z)   |
        \   |       2       \   |
         \  |                \  |
          \ |                 \ |
           \|                  \|
        (0,0,0) ----------- (x,0,0)

        """
        
        glBegin(GL_LINE_STRIP)
        glColor3f(0.0,0.0,1.0)
       
        # Face 1 
        glVertex3f(0,0,0)
        glVertex3f(0,0,z)
        glVertex3f(0,y,z)
        glVertex3f(0,y,0)
        glVertex3f(0,0,0)

        # Face 2
        glVertex3f(x,0,0)
        glVertex3f(x,y,0)
        glVertex3f(0,y,0)

        # Face 3
        glVertex3f(0,y,z)
        glVertex3f(x,y,z)
        glVertex3f(x,y,0)
        glVertex3f(0,y,0)

        # Face 5
        glVertex3f(0,y,z)
        glVertex3f(0,0,z)
        glVertex3f(x,0,z)
        glVertex3f(x,y,z)
        glVertex3f(0,y,z)

        # Face 6
        glVertex3f(x,y,z)
        glVertex3f(x,y,0)
        glVertex3f(x,0,0)
        glVertex3f(x,0,z)

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
            gluSphere(self.sphere,20.0,10,10)
            glPopMatrix()

            #glBegin(GL_POLYGON)
            #glColor3f(p.colour[0],p.colour[1],p.colour[2])
            #for angle in range(6):
            #    a = radians(angle*60)
            #    glVertex2f(r[0]+sin(a)*radius,r[1]+cos(a)*radius)
            #glEnd()

    def draw_neighbours(self,p):
        """ Issues the opengl commands to draw 
            lines connecting neighbours. """
        for i in range (p.nl_default.nip):
            glBegin(GL_LINES)
            glColor3f(0.3,0.3,0.3)
            a = p.nl_default.iap[i,0]
            b = p.nl_default.iap[i,1]
            glVertex3f(
                p.r[a,0]*self.xmap,
                p.r[a,1]*self.ymap,
                p.r[a,2]*self.zmap)
            glVertex3f(
                p.r[b,0]*self.xmap,
                p.r[b,1]*self.ymap,
                p.r[b,2]*self.zmap)
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
        glFlush()
        glPopMatrix()

        self.reset_perspective()

    def clear(self): 
        self.win.dispatch_events()
        glClear(GL_COLOR_BUFFER_BIT)
        glLoadIdentity()
        self.win.clear()
        
    def redraw(self,p):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        self.draw_neighbours(p)
        self.draw_particles(p)
        self.draw_box()
        self.hud(p)

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
        self.tbcam.update_modelview()

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

    #self.zoomin = 0.0
    #self.zoomout = 0.0
    #self.yawup = 0.0
    #self.yawdown = 0.0
    #self.panleft = 0.0
    #self.panright = 0.0

    #def update_eye(self):
    #        self.tbcam.cam_eye[2] += 100
    #        self.tbcam.update_modelview()

    def wsad_press(self,symbol,modifiers):
        """ zoom and rotate directly with the fps keys.
        """
        if symbol == pyglet.window.key.W:
            self.tbcam.cam_eye[2] += 100
            self.tbcam.update_modelview()
        if symbol == pyglet.window.key.S:
            self.tbcam.cam_eye[2] -= 100
            self.tbcam.update_modelview()
        if symbol == pyglet.window.key.A:
            self.tbcam.cam_eye[1] -= 30
        if symbol == pyglet.window.key.D:
            self.tbcam.cam_eye[1] += 30
        if symbol == pyglet.window.key.Q:
            self.tbcam.cam_eye[0] -= 30
        if symbol == pyglet.window.key.E:
            self.tbcam.cam_eye[0] += 30

    def wsad_release(self,symbol,modifiers):
        """ zoom and pan directly with the fps keys.
        """
        if symbol == pyglet.window.key.W:
            # zoom in
            self.tbcam.cam_eye[2] += 100
            self.tbcam.update_modelview()
        if symbol == pyglet.window.key.S:
            # zoom out
            self.tbcam.cam_eye[2] -= 100
            self.tbcam.update_modelview()
        if symbol == pyglet.window.key.A:
            self.tbcam.cam_eye[1] -= 30
            # pan left
            print 'you called pan left'
        if symbol == pyglet.window.key.D:
            self.tbcam.cam_eye[1] += 30
            # pan right
            print 'you called pan right'


class SmoothParticleView(ParticleView):
    """ Extends ParticleView to provide specialised visualisation of
        smooth particle systems. """
    def __init__(self):
        #config = Config(alpha_size=8)
        ParticleView.__init__(self)
        self.win = pyglet.window.Window(WINDOW_WIDTH,WINDOW_HEIGHT
            ,visible=False
            ,caption='Pyticles')
       
        self.win.set_visible()
       
        self.maxvol = pyglet.text.Label("maxvol label",font_name="Arial", \
            font_size=12,color =(255,0,0,255),x=200,y=18 )

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
        self.render_density(p)
        self.hud(p)
        self.draw_neighbours(p)
        self.draw_particles(p)
        #self.draw_particle_sprites(p)

