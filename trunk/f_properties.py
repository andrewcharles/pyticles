""" Wrapper around Fortran SPH code to compute SPH properties """
import sys
sys.path.append('/Users/acharles/masters/active/fsph')
import fkernel
import sphlib
from eos import vdw_co, vdw_hc
import numpy as np

def spam_properties(p,nl,h):
    """ Calculates and assigns:

        * kernel values
        * kernel gradient values
        * smoothed particle summation densities
        * velocity gradient

        Particle arrays will not be in fortran order.
        See if we can just do copies ($$$$$)

        nl.iap -- integer list of interactions
        nl.dv -- velocity difference

        p.rho -- particle density
        p.gradv -- particle velocity gradient
       
        does not take into account a longer cohesive smoothing length.
        this is a design problem for which I have notes already and
        will need to implement a neat solution. For now this is just to
        test the speed of raw python vs cython vs fortran

    """

    n = p.n
    d = p.dim
    ni = nl.nip

    v = np.reshape(p.v[0:n,:].copy(),(n,d),order='F')#.transpose()
    dv =  np.reshape(nl.dv[0:ni,:].copy(),(ni,d),order='F')#.transpose()
    nlist =  np.reshape(nl.iap[0:ni,:].copy(),(ni,2),order='F')#.transpose()
    
    #nlist = np.reshape(nlist,(ni,2),order='F')
    #nlist = nlist.transpose()
    #nlist = np.array(nlist[:,:].transpose(),order='F')

    T = p.t[0:n].reshape((n),order='F')
    T[:] = 1.0

    sml = np.reshape(p.h[0:n].copy(),(n,1),order='F')
    rij = np.reshape(nl.rij[0:ni].copy(),(ni,1),order='F')
    drij =np.reshape(nl.drij[0:ni,:].copy(),(ni,d),order='F')
    
    w =np.reshape( nl.wij[0:ni].copy(),(ni),order='F')
    dwdx = np.reshape(nl.dwij[0:ni,:].copy(),(ni,d),order='F')
    
    rho = np.reshape(p.rho[0:n].copy(),(n),order='F')
    gradv = np.reshape(p.gradv[0:n,:,:].copy(),(n,d,d),order='F')
    grad_rho = np.zeros((n,d),order='F')
    mass = np.reshape(p.m[0:n].copy(),(n,1),order='F')

    if ni == 0:
        # set values to isolated particle values
        # and issue a warning
        return

    #print nlist.flags

    # Velocity diff
    dv = sphlib.sphlib.calc_dv(dv,nlist+1,v)

    # Kernels and kernel gradients
    w,dwdx = fkernel.kernel.smoothing_kernels(rij,drij,nlist+1,sml,1)

    # Density summation
    fkernel.kernel.density_sum(rho,grad_rho,nlist,sml,mass,w,dwdx,1)

    # Grad v
    grad_v = np.zeros((n,d,d),order='F')
    sphlib.sphlib.calc_grad_v(grad_v,nlist,dwdx,dv,mass,rho)

    # (python implementation)
    phc = p.p[0:n]
    pco = p.pco[0:n]

    phc = vdw_hc(rho,T)
    pco = vdw_co(rho,T)

    # Resend data to python object
    p.rho[0:n] = rho[0:n]
    nl.wij[0:ni] = w[0:ni]
    nl.dwij[0:ni,:] = dwdx[0:ni,:]
    p.p[0:n] = phc[0:n]
    p.pco[0:n] = pco[0:n]
    






