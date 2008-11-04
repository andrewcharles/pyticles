"""
    Particle system model. Contains data structures and basic functions
    for a particle system and associated neighbour list.
    Copyright Andrew Charles 2008
    All rights reserved.
    This module is new BSD licensed.
"""

import random
import numpy
import math
import forces
import neighbour_list
import properties
import scipy
import integrator
import configuration
import box
#from Numeric import *
#import pdb

dt = 0.1
XMAX = 64 #64
YMAX = 32 #48
N = 25 
MAXN =100 
DIM = 2
VMAX = 0.1
CUTOFF_RADIUS =6 
VACUUM_VISCOSITY = 0.1
NVARS = 8 # this is now 8 with pco
RAMAL = 0.5  #Amalgamation radius
VSPLIT = 10.0 
rhosplit = 0.0
ADKE = True # Sigalotti style adaptive density kernel estimation
# make the mins just 0
AMALGAMATE = False
SPLIT = False

# variables for the integrator - put these somewhere cleaver
verbose = False


class Particles:
    """ A group of similar particles. 
        colour: a 3 tuple giving the RGB colour to render the particles in
    """
    def __init__(self,n):
        self.n = n
        self.maxn = MAXN
        self.dim = DIM 
        self.box = box.Box(self)
        self.V_split = VSPLIT
        self.r_amalg = RAMAL
        
        # basic mechanical properties
        self.m = numpy.zeros(self.maxn)
        self.m[:] = 1.
        self.h = numpy.zeros(self.maxn)
        self.h[:] = 3.
        self.v = VMAX * (numpy.random.random([self.maxn,2]) - 0.5)
        self.r = self.box.xmax * numpy.random.random([self.maxn,self.dim])
        self.rdot = numpy.zeros(self.r.shape)
        self.vdot = numpy.zeros(self.r.shape)
       
        #thermal properties
        self.t = numpy.ones(self.maxn)
        self.t[:] = 0.4
 
        #def grid(n,xside,yside,origin,spacing=1.0):
        self.r[0:self.n] = configuration.grid(self.n,5,5,(20,20),spacing=0.8)
        self.colour = 1.0,0.0,0.0 
        # sph properties
        self.rho = numpy.zeros(self.maxn)
        # pressure is just a scalar for now but this will
        # become the full pressure tensor in time.

        # I added another dimension to pressure, so that we have
        # the cohesive pressure and the repulsive pressure
        self.p = numpy.zeros([self.maxn,2])
        
        self.nlists = []
        # the next line just creates a placeholder list
        self.nl_default = neighbour_list.NeighbourList(self,10.0)
        
        self.rebuild_lists()
        for nl in self.nlists: 
            properties.spam_properties(self,nl,nl.cutoff_radius)

        self.x = numpy.zeros([NVARS,self.maxn])
        self.xdot = numpy.zeros([NVARS,self.maxn])


    def split(self,i):
        """ Splits particle i into a ragtag collection of four
            particles. Why four? Because calculating the new
            positions should be easy.

            The daughter particle is in the centre. The son particles
            are displaced slightly.

            See Feldman and Bonet, Dynamic refinement and boundary contact
            forces in SPH with applications in fluid flow problems.
            International Journal for Numerical Methods in Engineering.
            2007

        """
        alpha = 0.6
        eps = 2.6

        if self.n > self.maxn-3:
            print "cannot refine any further"
            return False
       
        # The son 
        self.m[i] = self.m[i] / 4.0
        #self.h[i] = self.h[i] * alpha

        # Daughter 1
        self.r[self.n] = self.r[i] + eps*numpy.array([0,1])
        self.m[self.n] = self.m[i] 
        self.v[self.n] = self.v[i]
        
        # Daughter 2
        self.r[self.n+1] = self.r[i] + eps*numpy.array([0.866025,-0.5])
        self.m[self.n+1] = self.m[i] 
        self.v[self.n+1] = self.v[i]
 
        # Daughter 3
        self.r[self.n+2] = self.r[i] + eps*numpy.array([-0.866025,-0.5])
        self.m[self.n+2] = self.m[i] 
        self.v[self.n+2] = self.v[i]
        
        self.n = self.n+3
        #print "There are now ",self.n,"particles"
        return True

    def create_particle(self,x,y):
        self.r[self.n] = x,y
        self.m[self.n] = self.m[self.n-1] 
        self.v[self.n] = 0,0
        self.n = self.n+1
        for nl in self.nlists: 
            nl.rebuild_list = True

    def check_refine(self):
        split = False
        for i in range(self.n):
            #V = self.m[i]/self.rho[i]
            #if V > self.V_split:
            if self.rho[i] < rhosplit:
                if verbose: print "V ",i," is ",V," - splitting"
                split=True
                self.split(i)        
            if split:
                for nl in self.nlists: 
                    nl.rebuild_list = True
              
    def amalgamate(self,i,j):
        """ Amalgamates particles i and j, merging them together to become
            one awesome robot with supernatural powers.
        """
        # conserve momentum
        self.v[i] = (self.v[i]*self.m[i]+self.v[j]*self.m[j])/ \
                    (self.m[i]+self.m[j])
        self.r[i] = (self.r[j] - self.r[i])/2 + self.r[j] 
        self.m[i] = self.m[i] + self.m[j]
        self.r[j] = self.r[self.n-1]
        self.v[j] = self.v[self.n-1]
        self.m[j] = self.m[self.n-1]
        self.n = self.n - 1

    def check_amalg(self,nl):
        for k in range(nl.nip):
            i = nl.iap[k,0]
            j = nl.iap[k,1]
            if nl.rij[k] < self.r_amalg:
                #print nl.rij[k]
                #print "amalgamating",i,j
                self.amalgamate(i,j)

    def rebuild_lists(self):
        """ rebuilds all nlists """
        
        for nl in self.nlists: 
            if nl.rebuild_list:
                nl.build_nl_verlet()

    def update(self):
        """ Update the particle system, using the
            neighbour list supplied.
        """
        # how to check that p is of class Particle?

        if SPLIT:
            self.check_refine()
        if AMALGAMATE:
            self.check_amalg(self.nl_default) 
        self.rebuild_lists()
        self.derivatives()
        
        for nl in self.nlists: 
            properties.spam_properties(self,nl,nl.cutoff_radius)
        # now integrate numerically
        integrator.rk4(self.gather_state,self.derivatives, \
                       self.gather_derivatives,self.scatter_state,dt)
        #integrator.euler(self.gather_state,self.derivatives, \
        #               self.gather_derivatives,self.scatter_state,dt)
        
        
        self.r[:,:] = self.r[:,:] + self.rdot[:,:]*dt
        self.v[:,:] = self.v[:,:] + self.vdot[:,:]*dt
        self.box.apply_mirror_bounds(self)
        #self.box.apply_periodic_bounds(self)
        
        #self.rebuild_lists()

    def gather_state(self):
        """ Maps the particle system to a state vector for integration
        """
        self.x[0,0:self.n] = self.m[0:self.n]
        self.x[1,0:self.n] = self.r[0:self.n,0]
        self.x[2,0:self.n] = self.r[0:self.n,1]
        self.x[3,0:self.n] = self.v[0:self.n,0]
        self.x[4,0:self.n] = self.v[0:self.n,1]
        self.x[5,0:self.n] = self.rho[0:self.n]
        self.x[6,0:self.n] = self.p[0:self.n,0]
        # added second component of pressure
        self.x[7,0:self.n] = self.p[0:self.n,1]
        return(self.x)

    def scatter_state(self,x):
        """ Maps the state vector to a particle system
        """
        self.m[0:self.n] = x[0,0:self.n] 
        self.r[0:self.n,0] = x[1,0:self.n]
        self.r[0:self.n,1] = x[2,0:self.n]
        self.v[0:self.n:,0] = x[3,0:self.n]
        self.v[0:self.n:,1] = x[4,0:self.n]
        self.rho[0:self.n] = x[5,0:self.n]
        self.p[0:self.n,0] = x[6,0:self.n]
        self.p[0:self.n,1] = x[7,0:self.n]

    def gather_derivatives(self):
        """ Maps particle system's derivatives to a state vector
        """
        self.xdot[0,0:self.n] = 0
        self.xdot[1,0:self.n] = self.v[0:self.n,0]
        self.xdot[2,0:self.n] = self.v[0:self.n,1]
        self.xdot[3,0:self.n] = self.vdot[0:self.n,0]
        self.xdot[4,0:self.n] = self.vdot[0:self.n,1]
        self.xdot[5,0:self.n] = 0
        self.xdot[6,0:self.n] = 0
        self.xdot[7,0:self.n] = 0
        return self.xdot

    def derivatives(self):
        """ get the rate of change of each variable 
            for every particle 
        """
        self.rdot = self.v
        # random velocity
        #self.vdot = numpy.random.random(self.v.shape)-0.5
        #hookes law
        self.vdot[:,:] = 0.0
    
        for nl in self.nlists: 
            if nl.rebuild_list:
                nl.build_nl_verlet()
            for force in nl.forces:
                force.apply()

   #     forces.apply_forces(self,nl)

