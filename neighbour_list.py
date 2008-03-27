"""
    Neighbour list class that works with particle and force
"""
import numpy


CUTOFF_RADIUS = 100
DIM = 2
XMAX = 20 #64
YMAX = 20 #48

class NeighbourList:
    """ A neighbour list that works with the Particle and Force class """
    """ nip: number of interacting paris
        max_interactions: self describing
        iap: current number of pairs
        drij: displacement between pairs
        rij: distance between pairs
        wij: interpolation kernel between pairs
        dwij: interpolation kernel gradient between pairs
    """
    def __init__(self,particle):
        self.nip = 0
        self.particle = particle
        self.cutoff_radius = CUTOFF_RADIUS
        self.cutoff_radius_sq = CUTOFF_RADIUS**2
        self.max_interactions = (particle.n * particle.n) / 2 - 1   
        self.iap = numpy.zeros((self.max_interactions,2),dtype=int)
        self.rij = numpy.zeros(self.max_interactions)
        self.drij = numpy.zeros((self.max_interactions,DIM))
        self.wij = numpy.zeros(self.max_interactions)
        self.dwij = numpy.zeros((self.max_interactions,DIM))

    def build_nl_brute_force(self):
        """ When in doubt, use brute force """
        i=0
        j=0
        self.nip = 0
        for i in range(self.particle.n):
            for j in range(i+1,self.particle.n):
               self.iap[self.nip,0] = i
               self.iap[self.nip,1] = j
               self.nip=self.nip+1

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
        k = 0
        for i in range(self.particle.n):
            for j in range(i+1,self.particle.n):
                self.drij[k,0] = self.particle.r[j,0] - self.particle.r[i,0]
                self.drij[k,1] = self.particle.r[j,1] - self.particle.r[i,1]
		self.minimum_image(self.drij[k,:],XMAX,YMAX)
                rsquared = self.drij[k,0]**2 + self.drij[k,1]**2
                if (rsquared < cutsq):
                    self.iap[k,0] = i
                    self.iap[k,1] = j
                    #self.nip=self.nip+1
                    k += 1
        self.nip = k     
            

