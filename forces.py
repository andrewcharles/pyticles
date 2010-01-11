""" 
    A generic force class ,and several force objects that operate on particle
    systems.

    Force -- A generic pairwise force that can be instantiated but
            does nothing.

    HookesForce -- a pairwise force that is linear in particle distance


    Andrew Charles
"""

import numpy
import math
import scipy
import neighbour_list
import sys
sys.path.append('/Users/acharles/masters/active/fsph')
from collision import collision


class Force:
    """ A generic pairwise particle force
        systems. The assumption is that the force is mutual. 
        p1: The first particle set
        nl: list of interacting pairs (p1,p2)

        This will become a subclass of Controller

    """ 

    def __init__(self,particles,nl,cutoff=100.0):
        self.p = particles
        self.nl = nl
        self.cutoff = cutoff
        self.cutoffsq = cutoff*cutoff
    
    def apply(self):
        """ Apply the force to all particles in the nlist """
        for k in range(self.nl.nip):
            if self.nl.rij[k]**2 <= self.cutoffsq:
                self.apply_force(k)

    def apply_sorted(self):
        """ Apply the force to all particles in the sorted nlist """
        """ Todo: check that the list is of type sorted!
        """
        for k in range(self.nl.nip):
            if self.nl.rij[k]**2 <= self.cutoffsq:
                self.apply_force(k)
            else:
                return


class HookesForce(Force):
    """
         Magnitude of force is k(|r-ro|)
         Assumes drij has been calculated and stored in nlist
         Takes a particle p and an nlist pair reference k
    """

    def __init__(self,particles,neighbour_list,cutoff=50.0,k=1.0,w=5.0):
        Force.__init__(self,particles,neighbour_list,cutoff=cutoff)
        self.w = w
        self.k = k

    def apply_force(self,k):
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p

        drx =  self.nl.drij[k,0] 
        dry =  self.nl.drij[k,1]
        rdist = self.nl.rij[k]
        rsquared = rdist**2 
        fmag = abs (self.k * ( rdist - self.w ) )
            
        #resolve into components
        dvx = fmag * ( drx ) / rdist
        dvy = fmag * ( dry ) / rdist
        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy


class CollisionForce(Force):
   
    def __init__(self,particles,neighbour_list,cutoff=5.0):
        Force.__init__(self,particles,neighbour_list,cutoff=cutoff)

    def apply_force(self,k):
        """ A hard collision is just an instantanous force. More
            like an impulse maybe. Anyhow, it changes v directly.
            TODO: rollback and roll forward to get rid of sticky balls.
        """
        debug = False
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        
        dr = self.nl.drij[k]
        drsq = dr[0]**2 + dr[1]**2
      
        # Divergence
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


class SpamForce2d(Force):
    
    def __init__(self,particles,neighbour_list,cutoff=5.0):
        Force.__init__(self,particles,neighbour_list,cutoff=cutoff)

    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        pri = self.p.p[i]
        prj = self.p.p[j]
        dwdx = self.nl.dwij[k,:]
       
        ps = (pri/p.rho[i]**2 + prj/p.rho[j]**2)

        dvx = ps * dwdx[0] 
        dvy = ps * dwdx[1] 

        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy

        du = 0.5 * (dvx * self.nl.dv[k,0] + dvy * self.nl.dv[k,1]) 
        p.udot[i] += du * p.m[j]
        p.udot[j] += du * p.m[i]


class CohesiveSpamForce2d(Force):
   
    def __init__(self,particles,neighbour_list,cutoff=10.0):
        Force.__init__(self,particles,neighbour_list,cutoff=cutoff)

    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """
        """p = self.p
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        pri = p.pco[i]
        prj = p.pco[j]
        dwdx = self.nl.dwij[k,:]
        dvx =  (pri/p.rho[i]**2 + prj/p.rho[j]**2) * dwdx[0]
        dvy =  (pri/p.rho[i]**2 + prj/p.rho[j]**2) * dwdx[1]
        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy"""

        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        pri = self.p.p[i]
        prj = self.p.p[j]
        dwdx = self.nl.dwij[k,:]
       
        ps = (pri/p.rho[i]**2 + prj/p.rho[j]**2)

        dvx = ps * dwdx[0] 
        dvy = ps * dwdx[1] 

        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy

        du = 0.5 * (dvx * self.nl.dv[k,0] + dvy * self.nl.dv[k,1]) 
        p.udot[i] += du * p.m[j]
        p.udot[j] += du * p.m[i]


class SpamForce(Force):
    
    def __init__(self,particles,neighbour_list,cutoff=5.0):
        Force.__init__(self,particles,neighbour_list,cutoff=cutoff)

    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        pri = self.p.p[i]
        prj = self.p.p[j]
        dwdx = self.nl.dwij[k,:]
        dv = self.nl.dv[k,:]
       
        ps = (pri/p.rho[i]**2 + prj/p.rho[j]**2)

        ax = ps * dwdx[0] 
        ay = ps * dwdx[1] 
        az = ps * dwdx[2] 

        p.vdot[i,0] += ax
        p.vdot[i,1] += ay
        p.vdot[i,2] += az
        p.vdot[j,0] -= ax
        p.vdot[j,1] -= ay
        p.vdot[j,2] -= az

        du = 0.5 * (ax * dv[0] + ay * dv[1] + az * dv[2]) 
        p.udot[i] += du * p.m[j]
        p.udot[j] += du * p.m[i]


class CohesiveSpamForce(Force):
    """ Moderately ridiculous that this is a seperate function, there is
        only one line difference..."""
   
    def __init__(self,particles,neighbour_list,cutoff=10.0):
        Force.__init__(self,particles,neighbour_list,cutoff=cutoff)

    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        pri = self.p.pco[i]
        prj = self.p.pco[j]
        dwdx = self.nl.dwij[k,:]
        dv = self.nl.dv[k,:]
       
        ps = (pri/p.rho[i]**2 + prj/p.rho[j]**2)

        ax = ps * dwdx[0] 
        ay = ps * dwdx[1] 
        az = ps * dwdx[2] 

        p.vdot[i,0] += ax
        p.vdot[i,1] += ay
        p.vdot[i,2] += az
        p.vdot[j,0] -= ax
        p.vdot[j,1] -= ay
        p.vdot[j,2] -= az

        du = 0.5 * (ax * dv[0] + ay * dv[1] + az * dv[2]) 
        p.udot[i] += du * p.m[j]
        p.udot[j] += du * p.m[i]





class Gravity(Force):
    """Controller that attracts other objects with an inverse square force.

       Acceleration of affected particles is computed as 

                dv/dt = (g*M)/r^2

       and directed towards the centre of the attractive domain.

       To do:
        - can we speed up the sqrt and vector ops?
        - look ahead one frame for position?

    """
    def __init__(self,particles,neighbour_list,cutoff=100.0,g=1.0):
        Force.__init__(self,particles,neighbour_list,cutoff=cutoff)
        self.g = g

    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]

        p = self.p

        drx =  self.nl.drij[k,0] 
        dry =  self.nl.drij[k,1]
        rdist = self.nl.rij[k]
        rsquared = rdist**2 

        fmag = self.g*p.m[i]*p.m[j]/rsquared

        #resolve into components
        dvx = fmag * ( drx ) / rdist
        dvy = fmag * ( dry ) / rdist
        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy
       

class FortranCollisionForce(Force):
    """ A wrapper for the fortran collide3d subroutine.
    """

    def __init__(self,particles,neighbour_list,cutoff=1.0):
        Force.__init__(self,particles,neighbour_list,cutoff=cutoff)
        self.nl = neighbour_list

    def apply_force(self,k):
        """ A hard collision is just an instantanous force.
        """
        nl = self.nl
        i = nl.iap[k,0]
        j = nl.iap[k,1]
        p = self.p
        collision.collide3d(p.v[i],p.v[j],p.m[i],p.m[j],nl.drij[k],nl.rsq[k])


