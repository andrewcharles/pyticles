""" 
    A generic force class ,and several force objects that operate on particle
    systems.
    Andrew Charles
"""

import numpy
import math
import scipy
import neighbour_list

# CONSTANTS
CUTOFF_RADIUS = 10
DIM = 2
COLLISION_RADIUS_SQ = 0.1
VACUUM_VISCOSITY = 0.1
spring_k = 10.1
rest_distance =  0.2

class Force:
    """ A generic pairwise particle force """ 

    # todo: add overloaded initialisation
    def __init__(self,particles,nl):
        self.p = particles
        self.nl = nl
    
    def apply(self):
        """ Apply the force to all particles in the nlist """
        for k in range(self.nl.nip):
            #i = self.nl.iap[k,0]
            #j = self.nl.iap[k,1]
            self.apply_force(k)

class HookesForce(Force):
    def __init__(self,particles,nl):
        self.p = particles
        self.nl = nl
        self.nl.cutoff_radius_sq = CUTOFF_RADIUS**2 

    def apply_force(self,k):
        """ Takes a particle p and an nlist pair reference k
        """
        # magnitude of force is k(|r-ro|)
        # should calculate this in nlist allocation
        # and use the stored value
        # actually, should calc it in a seperate method so
        # we can get new distances without getting a new
        # nlist completely
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        drx =  self.nl.drij[k,0] 
        dry =  self.nl.drij[k,1]
        rdist = self.nl.rij[k]
        rsquared = rdist**2 
        fmag = abs (spring_k * ( rdist - rest_distance ) )
            
        #resolve into components
        dvx = fmag * ( drx ) / rdist
        dvy = fmag * ( dry ) / rdist
        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy


class CollisionForce(Force):
    
    def __init__(self,particles,nl):
        self.p = particles
        self.nl = nl
        self.nl.build_nl_brute_force()
        self.nl.cutoff_radius_sq = COLLISION_RADIUS_SQ 

    def apply_force(self,k):
        """ A hard collision is just an instantanous force. More
            like an impulse maybe. Anyhow, it changes v directly.
            TODO: rollback and roll forward to get rid of sticky balls.
        """

        if self.nl.rij[k]**2 > self.nl.cutoff_radius_sq:
            print "toofar"
            return
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        debug = False 
        # Calc v dot r, magrhosq, and others
        dr = self.nl.drij[k]
        
        drsq = dr[0]**2 + dr[1]**2
        vidotr = p.v[i,0]*dr[0] + p.v[i,1]*dr[1]
        vjdotr = (p.v[j,0]*dr[0] + p.v[j,1]*dr[1])

        # If the particles are moving away from each other do nothing
        if (vidotr < 0) and (vjdotr >  0):
            return

        # Calculate tangential and normal components
        # of velocity
        vtix = (vidotr/drsq) * dr[0] 
        vtiy = (vidotr/drsq) * dr[1]
        vnix = p.v[i,0] - vtix
        vniy = p.v[i,1] - vtiy

        vtjx = (vjdotr/drsq) * dr[0] 
        vtjy = (vjdotr/drsq) * dr[1]
        vnjx = p.v[j,0] - vtjx
        vnjy = p.v[j,1] - vtjy

        if debug:
            vmagsqi = p.v[i,0]**2 + p.v[i,1]**2
            vmagsqj = p.v[j,0]**2 + p.v[j,1]**2
            vmagsqrot = (vtix+vnix)**2 + (vtiy+vniy)**2
            mom_i = math.sqrt(vmagsqi)
            mom_j = math.sqrt(vmagsqj)
            print "Before"
            print mom_i,mom_j,mom_i+mom_j
            #print vmagsqreg,vmagsqrot

        # Transfer tangential component of momentum
        tmp = vtix
        vtix = vtjx * (p.m[j]/p.m[i])
        vtjx = tmp * (p.m[i]/p.m[j])
        
        tmp = vtiy
        vtiy = vtjy * (p.m[j]/p.m[i])
        vtjy = tmp * (p.m[i]/p.m[j])

        # Convert back to xy frame
        p.v[i,0] = vtix + vnix 
        p.v[i,1] = vtiy + vniy
        p.v[j,0] = vtjx + vnjx
        p.v[j,1] = vtjy + vnjy
        
        if debug:
            vmagsqi = p.v[i,0]**2 + p.v[i,1]**2
            vmagsqj = p.v[j,0]**2 + p.v[j,1]**2
            mom_i = math.sqrt(vmagsqi)
            mom_j = math.sqrt(vmagsqj)
            print "After"
            print mom_i,mom_j,mom_i+mom_j
            exit()



class SpamForce(Force):
    
    def __init__(self,particles,nl):
        self.p = particles
        self.nl = nl
        self.nl.cutoff_radius_sq = CUTOFF_RADIUS**2

    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """
        #def spam(p,i,j,dwdx):
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        dwdx = self.nl.dwij[k,:]
        dvx =  -(p.p[i]/p.rho[i]**2 + p.p[j]/p.rho[j]**2) * dwdx[0]
        dvy =  -(p.p[i]/p.rho[i]**2 + p.p[j]/p.rho[j]**2) * dwdx[1]
        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy

        ### derived properties
        # force accumulator
        # internal energy
        # smoothing length
        # udot accumulator (sphforce(rho,q,p)
        # Q heat flux tensor
        # P pressure tensor ( eos(rho,t) )
