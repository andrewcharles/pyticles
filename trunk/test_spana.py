""" I have done without proper test scripts for too long. 

"""

import unittest

class NeighbourListTest(unittest.TestCase):
    """ 2. Test the neighbour lists. """
    """ The default neighbour list is a brute force list with no cutoff.
    """

    def test_nlist(self):
        import particles
        import neighbour_list
        n = 3
        p = particles.ParticleSystem(n,d=3,maxn=5)
        p.r[0,:] = (0.0,0.0,0.0) 
        p.r[1,:] = (1.0,0.0,0.0) 
        p.r[2,:] = (0.0,0.0,1.0) 
        nl = neighbour_list.NeighbourList(p)
        nl.build()
        nl.separations()
        self.assertEqual(nl.rij[0],1.0)

    def test_verlet_list(self):
        import particles
        import neighbour_list
        n = 3
        p = particles.ParticleSystem(n,d=3,maxn=5)
        p.r[0,:] = (0.0,0.0,0.0) 
        p.r[1,:] = (1.0,0.0,0.0) 
        p.r[2,:] = (0.0,0.0,1.0) 
        nl = neighbour_list.VerletList(p,cutoff=10,tolerance=2)
        nl.build()
        nl.compress()
        nl.separations()
        self.assertEqual(nl.rij[0],1.0)
        nl.ponder_rebuild()
        p.r[0,:] = (100.0,100.0,100.0)
        nl.ponder_rebuild()
        self.assertEqual(nl.rebuild_list,True)

class KernelTest(unittest.TestCase):
    """ 1. Test the smooth particle kernel functions. """

    def test_c_kernel(self):
        import c_test
        c_test.kernel_test()

    def test_py_kernel(self):
        import spkernel
        w, dwdx = spkernel.lucy_kernel(0.0,[0.0,0.0,0.0],2.0)
        print 'Python kernel at zero distance ',w


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
