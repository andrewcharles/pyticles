""" 
    A generic force class ,and several force objects that operate on particle
    systems.

    Most forces removed so I can concentrate on the Cython sematics

    The assumption is that two interacting particles are not from the
    same system.

    Andrew Charles
"""

import numpy as np
cimport numpy as np
import math
import scipy
import neighbour_list
import particles
#cimport _particles

# CONSTANTS
CUTOFF_RADIUS = 10
DIM = 2
COLLISION_RADIUS_SQ = 1.0
VACUUM_VISCOSITY = 0.1
spring_k = 1.1
rest_distance =  0.2

ctypedef np.float_t DTYPE_t
DTYPE = np.float

class Force:
    """ A generic pairwise particle force, between two particle
        systems. The assumption is that the force is mutual. 
        
        p1: The first particle set
        p2: The second particle set
        nl: list of interacting pairs (p1,p2)

        The design is wrong here:

    """ 

    def __init__(self,particles1,particles2,nl):
        self.p1 = particles1
        self.p2 = particles2
        self.nl = nl
    
    def apply(self):
        """ Apply the force to all particles in the nlist """
        for k in range(self.nl.nip):
            self.apply_force(k)


class SpamForce(Force):

    #cdef public _particles.SmoothParticleSystem p

    def __init__(self,particles,nl):
        self.p = particles
        self.nl = nl
        #self.nl.cutoff_radius_sq = CUTOFF_RADIUS**2

    def apply(self):
        """Iterate over the neighbour list and apply the force to all
        particles.
        """

        p = self.p
        nl = self.nl

        cdef np.ndarray[np.int_t,ndim=2,mode='c'] _iap 
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _dwdx
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _p
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _rho
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _vdot

        _iap = self.nl.iap.astype(np.int)
        _dwdx = nl.dwij.astype(np.float)
        _p = p.p.astype(np.float)
        _rho = p.rho.astype(np.float)
        _vdot = p.vdot.astype(np.float)

        cdef unsigned int nip = self.nl.nip
        cdef unsigned int i,j,k,n
        cdef float dvx, dvy, dvz

        n = p.n

        for k in xrange(nip):
            i = _iap[k,0]
            j = _iap[k,1]
            dvx =  (_p[i]/_rho[i]**2 + _p[j]/_rho[j]**2) * _dwdx[k,0]
            dvy =  (_p[i]/_rho[i]**2 + _p[j]/_rho[j]**2) * _dwdx[k,1]
            dvz =  (_p[i]/_rho[i]**2 + _p[j]/_rho[j]**2) * _dwdx[k,2]
            _vdot[i,0] += dvx
            _vdot[i,1] += dvy
            _vdot[i,2] += dvz
            _vdot[j,0] += -dvx
            _vdot[j,1] += -dvy
            _vdot[j,2] += -dvz

        p.rdot[0:n,:] = p.v[0:n,:]
        p.vdot[0:n,:] = _vdot[:,:]


    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """
        pass


class CohesiveSpamForce(Force):
    
    def __init__(self,particles,nl):
        self.p = particles
        self.nl = nl

    def apply(self):
        """Iterate over the neighbour list and apply the force to all
        particles.
        """

        p = self.p
        nl = self.nl

        cdef np.ndarray[np.int_t,ndim=2,mode='c'] _iap 
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _dwdx
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _p
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _rho
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _vdot

        _iap = self.nl.iap.astype(np.int)
        _dwdx = nl.dwij.astype(np.float)
        _p = p.pco.astype(np.float)
        _rho = p.rho.astype(np.float)
        _vdot = p.vdot.astype(np.float)

        cdef unsigned int nip = self.nl.nip
        cdef unsigned int i,j,k,n
        cdef float dvx, dvy, dvz
        n = p.n

        for k in xrange(nip):
            i = _iap[k,0]
            j = _iap[k,1]
            dvx =  (_p[i]/_rho[i]**2 + _p[j]/_rho[j]**2) * _dwdx[k,0]
            dvy =  (_p[i]/_rho[i]**2 + _p[j]/_rho[j]**2) * _dwdx[k,1]
            _vdot[i,0] += dvx
            _vdot[i,1] += dvy
            _vdot[j,0] += -dvx
            _vdot[j,1] += -dvy

        p.rdot[0:n,:] = p.v[0:n,:]
        p.vdot[0:n,:] = _vdot[:,:]


