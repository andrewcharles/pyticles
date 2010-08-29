"""
    Neighbour list class that works with particle and force
    Copyright Andrew Charles
    Written 2008
    All rights reserved.

    iap[] is always the list of particles for the purposes of interaction.
    todo: add PD's sophisticated minimum image for non-rectangular boxes?

"""
import numpy as np
from pairsep import pairsep

DIM = 3

class NeighbourList:
    """ A neighbour list that works with the Particle class.
        build() -- creates the pairs list and sets nip and rijsq
        nip -- number of interacting pairs
        max_interactions -- self describing
        iap[nip] -- pair indices
        rsq[nip]: squared distance between all list members
        drij[nip,d]: displacement between pairs
        dv[nip,d]: velocity difference between pairs
        rij[nip]: distance between pairs
        wij[nip]: interpolation kernel between pairs
        dwij[nip]: interpolation kernel gradient between pairs
    """
    def __init__(self,particle):
        self.nip = 0 
        self.particle = particle
        self.max_interactions = (particle.maxn * particle.maxn) / 2 - 1   
        self.iap = np.zeros((self.max_interactions,2),dtype=int)
        self.rij = np.zeros(self.max_interactions,dtype=float)
        self.rsq = np.zeros(self.max_interactions)
        self.drij = np.zeros((self.max_interactions,DIM),dtype=float)
        self.dv = np.zeros((self.max_interactions,DIM),dtype=float)
        self.wij = np.zeros(self.max_interactions)
        self.wij_lr = np.zeros(self.max_interactions)
        self.dwij = np.zeros((self.max_interactions,DIM))
        self.dwij_lr = np.zeros((self.max_interactions,DIM))
        self.rebuild_list = False
        self.nforce = 0
        self.forces = []

    def build(self):
        """ When in doubt, use brute force """
        i=0
        j=0
        self.nip = 0
        self.rebuild_list = False
        for i in range(self.particle.n):
            for j in range(i+1,self.particle.n):
                self.iap[self.nip,0] = i
                self.iap[self.nip,1] = j
                self.nip += 1

    def compress(self):
        """There is no concept of compression for a brute force list."""
        print 'Cannot compress a brute list'
        pass

    def separations(self):
        """ Computes the distance between pairs in the list and stores
            the result in the array rij, indexed by the same k that
            indexes the interacting pairs array iap.
        """
        for k in range(self.nip):
            i = self.iap[k,0]
            j = self.iap[k,1]
            self.drij[k,0] = self.particle.r[j,0] - self.particle.r[i,0]
            self.drij[k,1] = self.particle.r[j,1] - self.particle.r[i,1]
            self.drij[k,2] = self.particle.r[j,2] - self.particle.r[i,2]
            self.minimum_image(
                self.drij[k,:],self.particle.box.xmax,
                self.particle.box.ymax,
                self.particle.box.zmax)
            rsquared = self.drij[k,0]**2 + self.drij[k,1]**2 + self.drij[k,2]**2
            self.rij[k] = np.sqrt(rsquared)
            self.dv[k,0] = self.particle.v[j,0] - self.particle.v[i,0]
            self.dv[k,1] = self.particle.v[j,1] - self.particle.v[i,1]
            self.dv[k,2] = self.particle.v[j,2] - self.particle.v[i,2]
            #print i,j,self.particle.r[i,:],self.particle.r[j,:],self.drij[k,:]

    def find_pair(self,i,j):
        """ Finds the pair index k for the pair i,j if it exists.
            If it does not exist, returns -1.
            This is not optimised because the only use for it at the moment
            is in testing.
        """
        for k in range(self.nip):
            m = self.iap[k,0]
            n = self.iap[k,1]
            if( (i == m) and (j == n) ):
                return k
        return -1

    def apply_minimum_image(self):
        for k in range(self.nip):
            self.minimum_image(
                self.drij[k,:],self.particle.box.xmax,
                self.particle.box.ymax,
                self.particle.box.zmax)

    def minimum_image(self,dr,xmax,ymax,zmax):
        """ Applies the minimum image convention to the distance
            between all particles.
            # --- a permutation of dimensions --- #
        """
        drx,dry,drz = dr
        if (drx > xmax/2.):
            drx = drx - xmax
        if (dry > ymax/2.):
            dry = dry - ymax 	
        if (drz > zmax/2.):
            drz = drz - zmax
        if (drx < -xmax/2.):
            drx = drx + xmax
        if (dry < -ymax/2.):
            dry = dry + ymax 	
        if (drz < -zmax/2.):
            drz = drz + zmax
        dr[:] = (drx,dry,drz)
        #return dr


class VerletList(NeighbourList):
    """A brute force list with the added step of pruning pairs that
        are outside of the interaction range plus a tolerance.
        If pairs' separation changes by more than the tolerance
        the list needs to be rebuilt.

        We build a list in a number of steps. 

        First we build(), iterating over all possible pairs of 
        particles. If the squared distance is less than than the interaction
        cutoff plus the tolerance, the pair is added to the nlist.

        At any time we can compress(), which prunes additional pairs from
        the list.

        cutoff_radius_sq
        tolerance_sq
        dsqmax -- particle displacement threshold before list rebuild
        r_old -- particle positions at list build

        In the future I will
            add the concept of a neighbour list and an interaction list,
            but for now this simply computes the distances between pairs
            in the iap list.

    """
    def __init__(self,particle,cutoff=2.0,tolerance=1.0):
        NeighbourList.__init__(self,particle)
        self.cutoff_radius = cutoff 
        self.cutoff_radius_sq = cutoff**2
        self.tolerance_sq = tolerance * tolerance
        self.r_old = np.zeros(self.particle.r.shape)

    def build(self):
        """ Brute force with a cutoff radius. """
        i = 0
        j = 0
        k = 0
        self.nip = 0
        self.rebuild_list = False
        self.r_old[:,:] = self.particle.r[:,:]
        for i in range(self.particle.n):
            for j in range(i+1,self.particle.n):
                self.drij[k,0] = self.particle.r[j,0] - self.particle.r[i,0]
                self.drij[k,1] = self.particle.r[j,1] - self.particle.r[i,1]
                self.drij[k,2] = self.particle.r[j,2] - self.particle.r[i,2]
                self.minimum_image(
                    self.drij[k,:],self.particle.box.xmax,
                    self.particle.box.ymax,
                    self.particle.box.zmax)
                rsquared = self.drij[k,0]**2 + self.drij[k,1]**2 + self.drij[k,2]**2
                if (rsquared < (self.cutoff_radius_sq + self.tolerance_sq)):
                    self.drij[k,0] = self.drij[k,0]
                    self.drij[k,1] = self.drij[k,1]
                    self.drij[k,2] = self.drij[k,2]
                    self.dv[k,0] = self.particle.v[j,0] - self.particle.v[i,0]
                    self.dv[k,1] = self.particle.v[j,1] - self.particle.v[i,1]
                    self.dv[k,2] = self.particle.v[j,2] - self.particle.v[i,2]
                    self.rsq[k] = rsquared
                    self.iap[k,0] = i
                    self.iap[k,1] = j
                    k += 1
                    self.nip += 1

    def compress(self):
        """ Elimates pairs from the neigbour list that are outside
            the maximum interaction radius plus tolerance. 
            This is a separate function because we may want to compress several
            times before rebuilding the list.
            q -- new running total of pairs
        """
        q = 0
        #separations(self)
        for k in range(self.nip):
            i = self.iap[k,0]
            j = self.iap[k,1]
            self.drij[k,0] = self.particle.r[j,0] - self.particle.r[i,0]
            self.drij[k,1] = self.particle.r[j,1] - self.particle.r[i,1]
            self.drij[k,2] = self.particle.r[j,2] - self.particle.r[i,2]
            self.minimum_image(self.drij[k,:],self.particle.box.xmax,
                self.particle.box.ymax,
                self.particle.box.zmax)
            drx,dry,drz = self.drij[k,0],self.drij[k,1],self.drij[k,2]
            rsquared = drx**2 + dry**2 + drz**2
            if (rsquared < (self.cutoff_radius_sq + self.tolerance_sq)):
                self.drij[q,0] = drx
                self.drij[q,1] = dry
                self.drij[q,2] = drz
                self.dv[q,0] = self.dv[k,0]
                self.dv[q,1] = self.dv[k,1]
                self.dv[q,2] = self.dv[k,2]
                self.rsq[q] = rsquared
                self.iap[q,0] = i
                self.iap[q,1] = j
                q += 1
            self.nip = q
        self.ponder_rebuild()

    def ponder_rebuild(self):
        """ Ponder the decision of whether to rebuild the list.
            Set rebuild_list to True if any particle has moved
            further than the threshold distance
        """
        dr = self.r_old - self.particle.r
        rsquared = dr[:,0]**2 + dr[:,1]**2 + dr[:,2]**2
        dsq = np.max(rsquared)
        if dsq > self.tolerance_sq:
            self.rebuild_list = True
            
    def separations(self):
        """ Computes the distance between pairs in the list and stores
            the result in the array rij, indexed by the same k that
            indexes the interacting pairs array iap.
            #def pairsep(r,pairs,dr,ds):
            Redefined to just call the cython version!
        """
        pairsep(self)
        for k in range(self.nip):
            #i = self.iap[k,0]
            #j = self.iap[k,1]
            #drij[k,0] = self.drij[k,0]
            #drij[k,1] = self.drij[k,1]
            #drij[k,2] = self.drij[k,2]
            self.minimum_image(self.drij[k,:],self.particle.box.xmax,
                self.particle.box.ymax,
                self.particle.box.zmax)

class SmoothVerletList:
    """A Verlet list with smooth particle neighbourly properties,
        which are computed when the distances are computed.
    """

    def compute():
        pass

class SortedVerletList(VerletList):
    """ A sorted verlet neighbour list.
        idx is always an index.
    """

    def separations(self):
       pairsep(self)
       self.sort_by_r()

    def sort_by_r(self):
        """ Sorts the list by interparticle separation.
        """
        idx = np.argsort(self.rsq[:])
        idx = idx[::-1]
        self.rij = self.rij[idx] 
        self.iap = self.iap[idx,:]
        self.rij = self.rij[idx]
        self.rsq = self.rsq[idx]
        self.dv = self.dv[idx]
        self.drij = self.drij[idx,:]
        self.wij = self.wij[idx]
        self.dwij = self.dwij[idx,:]
        self.wij_lr = self.wij_lr[idx]
        self.dwij_lr = self.dwij_lr[idx,:]

    def sort_by_i(self):
        """ Sorted by i so we can only do say each particle's nearest
            couple of neighbours.
        """
        pass

    def build(self):
        """ Call the parent build, and then sort. """
        VerletList.build(self)
        self.sort_by_r()

    def compress(self):
        """ Call the parent compress, and then sort. """
        VerletList.compress(self)
        self.sort_by_r()


class CouplingList(NeighbourList):
    """ A neighbour list for interactions between two particle systems."""
    def __init__(self,particle_a,particle_b,cutoff):
        self.nip = 0
        self.particle_a = particle_a
        self.particle_b = particle_b
        self.cutoff_radius = cutoff 
        self.cutoff_radius_sq = cutoff**2
        # Compute brute force maximum interactions
        self.max_interactions = (particle_a.maxn * particle_b.maxn) / 2. - 1   
        self.iap = np.zeros((self.max_interactions,2),dtype=int)
        self.rij = np.zeros(self.max_interactions)
        self.drij = np.zeros((self.max_interactions,DIM))
        self.wij = np.zeros(self.max_interactions)
        self.dwij = np.zeros((self.max_interactions,DIM))
        self.rebuild_list = False
        self.nforce = 0
        self.forces = []

    def build(self):
        """ Not sure if verlet is the right term. We build a brute
            force list, and then eliminate pairs outside the
            interaction radius.
            Avoid function call overhead by just repeating the
            distance calculation code.
        """
        cutsq = self.cutoff_radius_sq
        self.rebuild_list = True 
        k = 0
        for i in range(self.particle.n):
            for j in range(self.particle2.n):
                drx = self.particle2.r[j,0] - self.particle.r[i,0]
                dry = self.particle2.r[j,1] - self.particle.r[i,1]
                drz = self.particle.r[j,2] - self.particle.r[i,2]
                rsquared = drx**2 + dry**2 + drz**2
                self.minimum_image(drij[k,:],self.particle.box.xmax,
                    self.particle.box.ymax,
                    self.particle.box.zmax) 
                
                if (rsquared < (cutsq + tolerance)):
                    self.drij[k,0] = drx
                    self.drij[k,1] = dry
                    self.drij[k,2] = drz
                    self.iap[k,0] = i
                    self.iap[k,1] = j
                    self.rij[k] = np.sqrt(rsquared)
                    k += 1
