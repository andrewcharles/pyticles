""" Objects for mutating particle systems. Pair potentials, body forces
    and even boundary conditions can all be modelled as controllers.


"""

import numpy

class Controller:
    """ Mutates one or more ParticleSystems.
    """

    def __init__(self):
        """ Initialise the controller with a reference to the
            particle system it is to mutate.

            p -- ParticleSystem object to be mutated.

        """
        self.groups = []
        return

    def bind_particles(self,p):
        self.groups.append(p)

    def apply(self):
        """All controllers must have an apply method."""
        print "Do nothing"
        return


class Force(Controller):
    """ This kind of controller mutates particle systems by modifying
        their rates of change. Generic Forces are abstract, we don't
        ever expect to create them unless for testing.
    """

    def apply(self):
        print "Do nothing force"
        return


class BodyForce(Force):
    """ A constant one body mechanical force, with a direction and a
        magnitude.
    """

    def __init__(self,dir=(0,1,0),mag=-1000):
        """
        p -- ParticleGroup
        dir -- direction, either a 2 or 3 tuple depending on the particle
                system's dimensionality

        """
        Controller.__init__(self)
        self.dir = dir
        self.mag = mag

    def bind_particles(self,p):
        self.groups.append(p)
        if p.dim == 2:
            self.dir = self.dir[0:1]

    def apply(self):
        """
        """
        for p in self.groups:
            for i in range(p.n):
                p.vdot[i,:] += (self.mag * numpy.array(self.dir))


class PairForce(Force):
    """ A force that operates between pairs of particles.
        This one needs a data structure to keep track of the
        pairs.
    """

    def __init__(self,nl):
        """ We only need to initialise with the neighbour list
            because the neighbour list refers to the particles
            it points to. There may be a clever way to implement
            this so that the pair force routines do not need to
            know if the two particles are from different systems
            or not.
        """
