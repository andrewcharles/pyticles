"""
    Particle system model. Contains data structures and basic functions
    for a particle system and associated neighbour list.
    Might be worth going back and abstracting some of the higher level stuff
    out, applying some good OO principles.
"""

import random
import numpy
import math
import forces
import neighbour_list
import properties
import scipy
import integrator
#from Numeric import *
#import pdb

dt = 0.001
XMAX = 20 #64
YMAX = 20 #48
N =2 
DIM = 2
VMAX = 0.1
CUTOFF_RADIUS = 100
COLLISION_RADIUS_SQ = 1
VACUUM_VISCOSITY = 0.1
NVARS = 7
# make the mins just 0

class Box:
    
    def __init__(self,p):
        self.p = p

    def min_image(r):
        print "does not exist yet"

    def apply_periodic_bounds(self,p):
         """ applies periodic boundaries """
         for i in range(p.n):
            if p.r[i,0] > XMAX:
                p.r[i,0] = 0
            if p.r[i,0] < 0:
                p.r[i,0] = XMAX
            if p.r[i,1] > YMAX:
                p.r[i,1] = 0
            if p.r[i,1] < 0:
                p.r[i,1] = YMAX

    def apply_mirror_bounds(self,p):
         """ applies periodic boundaries """
         for i in range(p.n):
            if p.r[i,0] > XMAX:
                p.r[i,0] = XMAX
                p.v[i,0] = -p.v[i,0]
            if p.r[i,0] < 0:
                p.r[i,0] = 0
                p.v[i,0] = -p.v[i,0]
            if p.r[i,1] > YMAX:
                p.r[i,1] = YMAX 
                p.v[i,1] = -p.v[i,1]
            if p.r[i,1] < 0:
                p.r[i,1] = 0
                p.v[i,1] = -p.v[i,1]


class Particles:
    """ A group of similar particles """
    def __init__(self):
        self.n = N
        self.dim = DIM 
        self.box = Box(self)
        
        # basic mechanical properties
        self.m = numpy.zeros(self.n)
        self.m[:] = 1
        self.v = VMAX * (numpy.random.random([self.n,2]) - 0.5)
        self.r = XMAX * numpy.random.random([self.n,self.dim])
        self.rdot = numpy.zeros(self.r.shape)
        self.vdot = numpy.zeros(self.r.shape)
        
        # sph properties
        self.rho = numpy.zeros(self.n)
        # pressure is just a scalar for now but this will
        # become the full pressure tensor in time.
        self.p = numpy.zeros(self.n)
        
        self.nforce = 0
        self.forces = []

    def add_force(self,f):
        """ Adds a force to a particle system """
        self.forces.append(f)
        self.nforce += 1

    def update(self,nl):
        """ Update the particle system, using the
            neighbour list supplied.
        """
        # how to check that p is of class Particle?
        properties.spam_properties(self,nl)
        self.derivatives()
        
        # now integrate numerically
        integrator.imp_euler(self.gather_state,self.derivatives,self.gather_derivatives,self.scatter_state,dt)
        
    #    self.r[:,:] = integrator.euler()
    #    self.v[:,:] = integrator.euler()
        
        self.r[:,:] = self.r[:,:] + self.rdot[:,:]*dt
        self.v[:,:] = self.v[:,:] + self.vdot[:,:]*dt
        #self.box.apply_mirror_bounds(self)
        self.box.apply_periodic_bounds(self)

    def gather_state(self):
        """ Maps the particle system to a state vector for integration
        """
        x = numpy.zeros([NVARS,self.n])
        x[0,:] = self.m[:]
        x[1,:] = self.r[:,0]
        x[2,:] = self.r[:,1]
        x[3,:] = self.v[:,0]
        x[4,:] = self.v[:,1]
        x[5,:] = self.rho[:]
        x[6,:] = self.p[:]
        return(x)

    def scatter_state(self,x):
        """ Maps the state vector to a particle system
        """
        self.m = x[0,:] 
        self.r[:,0] = x[1,:]
        self.r[:,1] = x[2,:]
        self.v[:,0] = x[3,:]
        self.v[:,1] = x[4,:]
        self.rho = x[5,:]
        self.p = x[6,:]

    def gather_derivatives(self):
        """ Maps particle system's derivatives to a state vector
        """
        xdot = numpy.zeros([NVARS,self.n])
        xdot[0,:] = 0
        xdot[1,:] = self.v[:,0]
        xdot[2,:] = self.v[:,1]
        xdot[3,:] = self.vdot[:,0]
        xdot[4,:] = self.vdot[:,1]
        xdot[5,:] = 0
        xdot[6,:] = 0
        return xdot

    def derivatives(self):
        """ get the rate of change of each variable 
            for every particle 
        """
        self.rdot = self.v
        # random velocity
        #self.vdot = numpy.random.random(self.v.shape)-0.5
        #hookes law
        self.vdot[:,:] = 0.0
    
        for force in self.forces:
            # check rebuild condition and rebuild if needed
            force.nl.build_nl_verlet()
            force.apply()
   #     forces.apply_forces(self,nl)

