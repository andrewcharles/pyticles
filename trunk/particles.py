"""
    Particle system model. Contains data structures and basic functions
    for a particle system and associated neighbour list.

    Designed so that the neighbour list and force
    modules are accessed through here.

    I have a particle system class and a smooth particle system
    class because I am interested in the bits that are needed
    for a minimal particle system, and the bits that are added to
    make it an sph system.

    Copyright Andrew Charles 2008
    All rights reserved.
    This module is new BSD licensed.
"""

import random
import numpy as np
import math
#import forces
import neighbour_list
#import properties

# I don't want to have to select this here!
#import c_properties as properties
import f_properties as properties
from properties import vdw_energy,spam_properties,hamiltonian
import scipy
import integrator
import configuration
import box
import sys
#from Numeric import *
#import pdb
from integrator import rk4, euler, imp_euler
from time import time

#dt = 0.05

# variables for the integrator - put these somewhere cleaver
verbose = False

XMAX = 5 #64
YMAX = 5 #48
ZMAX = 1
VMAX = 0.1
CUTOFF_RADIUS = 6 
VACUUM_VISCOSITY = 0.1
RAMAL = 0.5  #Amalgamation radius
VSPLIT = 10.0 
rhosplit = 0.0
ADKE = True # Sigalotti style adaptive density kernel estimation
AMALGAMATE = False
SPLIT = False
ADVECTIVE = False


class ParticleSystem:
    """ A group of similar particles with basic mechanical properties.

    """
    def __init__(self,n,
            d=3,
            maxn=125,
            controllers=[],
            xmax=XMAX,
            ymax=YMAX,
            zmax=ZMAX,
            vmax=VMAX,
            simbox=None,
            rinit=None,
            side=(5,5,5),
            spacing=0.1):
        """
        DIMENSIONS
        n -- initial number of particles.
        maxn -- maximum number of particles.
        dim -- number of spatial dimensions (1-3).
        steps -- number of steps taken
        nlists -- neighbour lists associated with this particle system.
        colour -- a 3 tuple giving the particles' RBG color.
        box -- the simulation box. Should replace this with constraint forces.
        nlists -- neighbour lists associated with this system.
        forces -- internal forces associated with this system.
        configuration -- initial positions.
        side -- side lengths (only for grid positions)

        """
        self.n = n
        self.dim = d
        self.maxn = maxn
        self.dt = 0.0
        self.steps = 0

        # Basic mechanical properties
        if not simbox:
            self.box = box.MirrorBox(self,xmax=xmax,ymax=ymax,zmax=zmax)
        
        # Start with a random configuration
        self.r = self.box.xmax * np.random.random([self.maxn,self.dim])

        if rinit == 'grid':
           self.r[0:n,:] = configuration.grid3d(n,side,(xmax/2.,ymax/2.,zmax/2.)
                ,spacing=spacing)
    
        self.m = np.zeros(self.maxn)
        self.v = vmax * (np.random.random([self.maxn,self.dim]) - 0.5)
        self.rdot = np.zeros(self.r.shape)
        self.vdot = np.zeros(self.v.shape)
        self.mdot = np.zeros(self.m.shape)

        # Initialise values
        #self.r[0:self.n]=configuration.grid3d(self.n,5,5,(20,20,20),spacing=0.8)
        self.m[:] = 1.
        self.colour = 1.0,0.0,0.0 

        # State vectors to pass to numerical integrators
        n_variables = 7
        self.x = np.zeros([n_variables,self.maxn])
        self.xdot = np.zeros([n_variables,self.maxn])

        self.nlists = []
        self.forces = []

        self.controllers = controllers
        for controller in self.controllers:
            controller.bind_particles(self)

        """ Variables for measuring performance. """
        self.timing = {}
        self.timing['force time'] = -1
        self.timing['deriv time'] = -1
        self.timing['pairsep time'] = -1
        self.timing['update time'] = -1
        self.timing['integrate time'] = -1

    def create_particle(self,r,v=(0.0,0.0,0.0)):
        """Adds a new particle to the system.
        """
        self.r[self.n] = r
        self.m[self.n] = self.m[self.n-1] 
        self.v[self.n] = v
        self.n = self.n+1
        self.rebuild_lists()

    def rebuild_lists(self):
        """ rebuilds all nlists """
        
        for nl in self.nlists: 
            if nl.rebuild_list:
                nl.build()

    def update(self,dt):
        """ Update the particle system, using the
            neighbour list supplied.
        """
        self.rebuild_lists()
        imp_euler(self.gather_state,self.derivatives, \
                       self.gather_derivatives,self.scatter_state,dt)
        self.box.apply(self)
        self.steps += 1

# Right now these are hard coded to 3d. I am still mulling over the
# best approach the being able to do 2 and 1d if I want to.

    def gather_state(self):
        """ Maps the particle system to a state vector for integration
        """
        self.x[0,0:self.n] = self.m[0:self.n]
        self.x[1,0:self.n] = self.r[0:self.n,0]
        self.x[2,0:self.n] = self.r[0:self.n,1]
        self.x[3,0:self.n] = self.r[0:self.n,2]
        self.x[4,0:self.n] = self.v[0:self.n,0]
        self.x[5,0:self.n] = self.v[0:self.n,1]
        self.x[6,0:self.n] = self.v[0:self.n,2]
        return(self.x)

    def scatter_state(self,x):
        """ Maps the state vector to a particle system
        """
        self.m[0:self.n] = x[0,0:self.n] 
        self.r[0:self.n,0] = x[1,0:self.n]
        self.r[0:self.n,1] = x[2,0:self.n]
        self.r[0:self.n,2] = x[3,0:self.n]
        self.v[0:self.n:,0] = x[4,0:self.n]
        self.v[0:self.n:,1] = x[5,0:self.n]
        self.v[0:self.n:,2] = x[6,0:self.n]

    def gather_derivatives(self):
        """ Maps particle system's derivatives to a state vector
        """
        self.xdot[0,0:self.n] = self.mdot[0:self.n] 
        self.xdot[1,0:self.n] = self.rdot[0:self.n,0]
        self.xdot[2,0:self.n] = self.rdot[0:self.n,1]
        self.xdot[3,0:self.n] = self.rdot[0:self.n,2]
        self.xdot[4,0:self.n] = self.vdot[0:self.n,0]
        self.xdot[5,0:self.n] = self.vdot[0:self.n,1]
        self.xdot[6,0:self.n] = self.vdot[0:self.n,2]
        return self.xdot

    def derivatives(self):
        """ Compute the rate of change of each variable 
            for every particle. The force.apply() call
            accumulates the forces.
        """
        self.rdot = self.v
        self.vdot[:,:] = 0.0
    
        for nl in self.nlists: 
            nl.separations()
        
        for force in self.forces:
            force.apply()

        # Controllers is the new implementation of forces
        for controller in self.controllers:
            controller.apply()


class SmoothParticleSystem(ParticleSystem):
    """A particle system with additional properties to solve smooth
    particle equations of motion.

    This is actually a long-range/short range smooth particle system.
    I'd like to make the way this is structured more elegant - would
    like to have a simple smooth particle system, and have the long-short
    one inherit from it. Another solution is to have a smooth properties
    class and let the smooth particle system have multiple instances.

    As I don't know the best solution I am going to think on it.
    """

    def __init__(self,n,
            d=3,
            maxn=100,
            controllers=[],
            xmax=XMAX,
            ymax=YMAX,
            zmax=ZMAX,
            vmax=VMAX,
            rinit=None,
            side=None,
            spacing=None,
            temperature=1.0,
            hshort=1.0,
            hlong=3.0,
            simbox=None):
        ParticleSystem.__init__(self,n=n,d=d,
            xmax=xmax,
            ymax=ymax,
            zmax=zmax,
            rinit=rinit,
            side=side,
            spacing=spacing,
            vmax=vmax,
            maxn=maxn)
        self.V_split = VSPLIT
        self.r_amalg = RAMAL

        """
        SPH Properties
        --------------
        rho -- mass density
        rhodot -- time rate of change of mass density
        gradv --  spatial gradient of velocity
        t -- temperature
        h -- smoothing length
        p -- isotropic pressure (repulsive)
        pco -- isotropic pressure (cohesive)
        jq -- heat flux

        P -- pressure tensor

        timing -- a dictionary of average execution times for
                  particular subroutines.

        """

        self.rho = np.zeros(self.maxn)
        self.rhodot = np.zeros(self.rho.shape)
        self.gradv = np.zeros([self.maxn,self.dim,self.dim])
        self.jq = np.zeros([self.maxn,self.dim])
        #thermal properties
        self.t = np.ones(self.maxn)
        self.t[:] = temperature

        self.u = np.ones(self.maxn,dtype=float)
        self.udot = np.zeros(self.maxn,dtype=float)

        self.h = np.zeros(self.maxn)
        self.hlr = np.zeros(self.maxn)
        self.h[:] = hshort
        self.hlr[:] = hlong
        self.p = np.zeros([self.maxn])
        self.pco = np.zeros([self.maxn])
        self.P = np.zeros([self.maxn,self.dim,self.dim])
        
        #for nl in self.nlists: 
        #    nl.build()
        #    properties.spam_properties(self,nl,nl.cutoff_radius)

        #properties.spam_properties(self,self.nlists[0],self.nlists[1] \
        #    ,self.h[0:self.n],self.hlr[0:self.n])

        n_variables = 11
        self.x = np.zeros([n_variables,self.maxn])
        self.xdot = np.zeros([n_variables,self.maxn])
        self.timing['SPAM time'] = -1


    def split(self,i):
        """ Splits particle i into a ragtag collection of four
            particles. Why four? Because calculating the new
            positions is straightforward.

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
        self.r[self.n] = self.r[i] + eps*np.array([0,1])
        self.m[self.n] = self.m[i] 
        self.v[self.n] = self.v[i]
        
        # Daughter 2
        self.r[self.n+1] = self.r[i] + eps*np.array([0.866025,-0.5])
        self.m[self.n+1] = self.m[i] 
        self.v[self.n+1] = self.v[i]
 
        # Daughter 3
        self.r[self.n+2] = self.r[i] + eps*np.array([-0.866025,-0.5])
        self.m[self.n+2] = self.m[i] 
        self.v[self.n+2] = self.v[i]
        
        self.n = self.n+3
        #print "There are now ",self.n,"particles"
        return True


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
                nl.build()

    def update(self,dt):
        """ Update the particle system, using the
            neighbour list supplied.
        """
        t1 = time()

        if SPLIT:
            self.check_refine()
        if AMALGAMATE:
            self.check_amalg(self.nl_default)

        t = time()
        self.rebuild_lists()
        self.timing['nlist rebuild time'] = time() - t
        
        t = time()
        self.derivatives()
        self.timing['deriv time'] = time() - t
       
        t = time()
        imp_euler(self.gather_state,self.derivatives, \
            self.gather_derivatives,self.scatter_state,dt)
        self.timing['integrate time'] = time() - t
        
        self.box.apply(self)
        
        self.timing['update time'] = time() - t1
        self.steps += 1

    def gather_state(self):
        """ Maps the particle system to a state vector for integration
        """
        self.x[0,0:self.n] = self.m[0:self.n]
        self.x[1,0:self.n] = self.r[0:self.n,0]
        self.x[2,0:self.n] = self.r[0:self.n,1]
        self.x[3,0:self.n] = self.r[0:self.n,2]
        self.x[4,0:self.n] = self.v[0:self.n,0]
        self.x[5,0:self.n] = self.v[0:self.n,1]
        self.x[6,0:self.n] = self.v[0:self.n,2]
        self.x[7,0:self.n] = self.rho[0:self.n]
        self.x[8,0:self.n] = self.p[0:self.n]
        # added second component of pressure
        self.x[9,0:self.n] = self.pco[0:self.n]
        self.x[10,0:self.n] = self.u[0:self.n]
        return(self.x)

    def scatter_state(self,x):
        """ Maps the state vector to a particle system
        """
        self.m[0:self.n] = x[0,0:self.n] 
        self.r[0:self.n,0] = x[1,0:self.n]
        self.r[0:self.n,1] = x[2,0:self.n]
        self.r[0:self.n,2] = x[3,0:self.n]
        self.v[0:self.n:,0] = x[4,0:self.n]
        self.v[0:self.n:,1] = x[5,0:self.n]
        self.v[0:self.n:,2] = x[6,0:self.n]
        self.rho[0:self.n] = x[7,0:self.n]
        self.p[0:self.n] = x[8,0:self.n]
        self.pco[0:self.n] = x[9,0:self.n]
        self.u[0:self.n] = x[10,0:self.n]

    def gather_derivatives(self):
        """ Maps particle system's derivatives to a state vector
        """
        self.xdot[0,0:self.n] = self.mdot[0:self.n] 
        self.xdot[1,0:self.n] = self.rdot[0:self.n,0]
        self.xdot[2,0:self.n] = self.rdot[0:self.n,1]
        self.xdot[3,0:self.n] = self.rdot[0:self.n,2]
        self.xdot[4,0:self.n] = self.vdot[0:self.n,0]
        self.xdot[5,0:self.n] = self.vdot[0:self.n,1]
        self.xdot[6,0:self.n] = self.vdot[0:self.n,2]
        self.xdot[7,0:self.n] = self.rhodot[0:self.n] 
        self.xdot[8,0:self.n] = 0
        self.xdot[9,0:self.n] = 0
        self.xdot[10,0:self.n] = self.udot[0:self.n]
        return self.xdot

    def derivatives(self):
        """ get the rate of change of each variable 
            for every particle 
        """
        self.rdot = self.v
        # random velocity
        #self.vdot = np.random.random(self.v.shape)-0.5
        self.vdot[:,:] = 0.0
        self.udot[:] = 0.0

        t = time()
        for nl in self.nlists: 
            nl.separations()
    
        self.timing['pairsep time'] = (time() - t)

        t = time()

        properties.spam_properties(self,self.nlists[0] \
            ,self.h[0:self.n],self.hlr[0:self.n])
        
        self.timing['SPAM time'] = time() - t
        
        t = time()
        for force in self.forces:
            force.apply()
        self.timing['force time'] = time() - t
        
        if ADVECTIVE:
            self.rdot[:,:] = 0.0

#        print self.p[0],self.pco[0],self.t[0],self.u[0]
#        hamiltonian(self)
#        print np.mean(self.t)

