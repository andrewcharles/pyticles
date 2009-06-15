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
        print 'Number of pairs for',n,'particles',nl.nip
        print 'Particle positions',p.r[0:n]
        print 'Pair separations',nl.rij[0:nl.nip]
        self.assertEqual(nl.rij[0],1.0)

class KernelTest(unittest.TestCase):
    """ 1. Test the smooth particle kernel functions. """

    def test_c_kernel(self):
        import c_test
        c_test.kernel_test()

    def test_py_kernel(self):
        import spkernel
        w, dwdx = spkernel.lucy_kernel(0.0,[0.0,0.0,0.0],2.0)
        print 'Python kernel at zero distance ',w




if __name__=='__main__':
	unittest.main()
