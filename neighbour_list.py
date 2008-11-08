"""
    Neighbour list class that works with particle and force
    Copyright Andrew Charles 2008
    All rights reserved.
"""
import numpy


DIM = 2

class NeighbourList:
    """ A neighbour list that works with the Particle and Force class """
    """ nip: number of interacting pairs
        max_interactions: self describing
        iap: current number of pairs
        drij: displacement between pairs
        rij: distance between pairs
        wij: interpolation kernel between pairs
        dwij: interpolation kernel gradient between pairs
    """
    def __init__(self,particle,cutoff,particle2=None):
        self.nip = 0
        self.particle = particle
        if particle2 is not None:
            self.particle2 = particle2
        else:
            self.particle2 = None
        self.cutoff_radius = cutoff 
        self.cutoff_radius_sq = cutoff**2
        
        if particle2 is not None:
            self.max_interactions = (particle.maxn * particle2.maxn) / 2 - 1   
        self.max_interactions = (particle.maxn * particle.maxn) / 2 - 1   
        
        self.iap = numpy.zeros((self.max_interactions,2),dtype=int)
        self.rij = numpy.zeros(self.max_interactions)
        self.drij = numpy.zeros((self.max_interactions,DIM))
        self.wij = numpy.zeros(self.max_interactions)
        self.dwij = numpy.zeros((self.max_interactions,DIM))
        self.rebuild_list = False
        self.nforce = 0
        self.forces = []

    def build_nl_brute_force(self):
        """ When in doubt, use brute force """
        i=0
        j=0
        self.nip = 0
        k = 0
        self.rebuild_list = False
        for i in range(self.particle.n):
            for j in range(i+1,self.particle.n):
                self.drij[k,0] = self.particle.r[j,0] - self.particle.r[i,0]
                self.drij[k,1] = self.particle.r[j,1] - self.particle.r[i,1]
                rsquared = self.drij[k,0]**2 + self.drij[k,1]**2
                self.rij[k] = numpy.sqrt(rsquared)
                self.iap[self.nip,0] = i
                self.iap[self.nip,1] = j
                self.nip=self.nip+1
                k = self.nip

    def add_force(self,f):
        """ Adds a force to a neighbour list """
        self.forces.append(f)
        self.nforce += 1
    
    def minimum_image(self,dr,xmax,ymax):
        """ applies the minimum image convention to the distance
            between two particles. Based on code by Peter Daivis
            <insert reference to paper>
        """
        # don't use the general shape pbs
        # just the 2d rectangle
        if (dr[0] > xmax):
            dr[0] = dr[0] - xmax
        if (dr[1] > ymax):
            dr[1] = dr[1] - ymax 	

    def build_nl_verlet(self):
        """ Not sure if verlet is the right term. We build a brute
            force list, and then eliminate pairs outside the
            interaction radius.
            Avoid function call overhead by just repeating the
            distance calculation code.
        """
        cutsq = self.cutoff_radius_sq
        self.rebuild_list = True 
        k = 0
        if self.particle2 is None:
            for i in range(self.particle.n):
                for j in range(i+1,self.particle.n):
                    self.drij[k,0] = self.particle.r[j,0] - self.particle.r[i,0]
                    self.drij[k,1] = self.particle.r[j,1] - self.particle.r[i,1]
                    #self.minimum_image(self.drij[k,:],XMAX/2,YMAX/2)
                    rsquared = self.drij[k,0]**2 + self.drij[k,1]**2
                    if (rsquared < cutsq):
                        self.iap[k,0] = i
                        self.iap[k,1] = j
                        self.rij[k] = numpy.sqrt(rsquared)
                        k += 1

        else:
            for i in range(self.particle.n):
                for j in range(self.particle2.n):
                    self.drij[k,0] = self.particle2.r[j,0] - self.particle.r[i,0]
                    self.drij[k,1] = self.particle2.r[j,1] - self.particle.r[i,1]
                    #self.minimum_image(self.drij[k,:],XMAX/2,YMAX/2)
                    rsquared = self.drij[k,0]**2 + self.drij[k,1]**2
                    if (rsquared < cutsq):
                        self.iap[k,0] = i
                        self.iap[k,1] = j
                        self.rij[k] = numpy.sqrt(rsquared)
                        k += 1



        self.nip = k     
            


