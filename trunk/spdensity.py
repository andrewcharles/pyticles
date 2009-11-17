""" A set of functions for computing smoothed particle summation
    densities.
"""
import sys
sys.path.append("../pyticles")
import spkernel
import numpy as np

def sum_pdensity(r,m):
    """ This computes summation density.
        pairs is an array giving the i,j
        indices of interacting pairs
        r is the particle positions 
        m is the particle masses
        We make a neighbour list and calculate the kernel
        in house.
    """
    #for i in range(n):
    print 'wtf'
    #    nlist = neighbours.Nlist(n,n)
    #    nlist.build(r[i],r,h,periodic=True,xmax=xmax,ymax=ymax)
    #    w = numpy.zeros(len(nlist.nrefs))
    #    for j in range(len(nlist.nrefs)): 
    #        w[j] = spkernel.kernel(nlist.ds[j],nlist.dr[j],h[j],'lucy')[0]
    #    rho[i] = spdensity.summation_density(w,m)
    #return rho

def summation_density(w,m):
    """ 
        w - a numpy array of kernel values
        m - a numpy array of masses
        function returns the value of the density
    """
    n = w.size
    rho = 0.0
    for i in range(n):
        rho += m[i] * w[i]
    return rho 

def sum_density(m,h,nl,rij,drij):
    zdist = np.zeros(drij[0,:].shape)
    zerokern = spkernel.lucy_kernel(0.0,zdist,h[0])[0]
    n = h.size
    ni = nl.shape[0]
    rho = np.zeros(n)
    rho[:] = zerokern

    # calc the kernels, velocities and densities
    for k in range(ni):
        i = nl[k,0]
        j = nl[k,1]
        wij = spkernel.lucy_kernel(rij[k],drij[k,:],h[i])[0]
        rho[i] += wij * m[j]
        rho[j] += wij * m[i]

    return rho
