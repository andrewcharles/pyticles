import unittest

class ForceTest(unittest.TestCase):
    """ Test the sph force subroutine. """

    def test_spam_force(self):
        import particles
        import neighbour_list
        import properties
        import forces
        n = 2
        p = particles.SmoothParticleSystem(n,d=3,maxn=5)
        p.r[0,:] = (0.0,0.0,0.0) 
        p.r[1,:] = (1.0,0.0,0.0) 
        nl = neighbour_list.VerletList(p,cutoff=10,tolerance=2)
        nl.build()
        nl.compress()
        nl.separations()
        properties.spam_properties(p,nl,5.0)
        f = forces.SpamForce(p,nl)
        f.apply()


if __name__=='__main__':
	unittest.main()
