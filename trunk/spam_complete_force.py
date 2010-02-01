from forces import Force

""" A wrapper around the fortran sphforce3d routine.
    Computes some sph properties first.

"""

import sys
sys.path.append('/Users/acharles/masters/active/fsph')
import fkernel
import sphlib
import feos
import numpy as np
import sphforce3d
from time import time

feos.eos.adash = 2.0
feos.eos.bdash = 0.5
feos.eos.kbdash = 1.0
sphforce3d.sigma = 1.0
sphforce3d.rcoef = 0.5
sphforce3d.cgrad = 0.0

class SpamComplete(Force):
    """ Compute all the smooth particle properties and forces
        in a single function that uses the fortran routines.
    """

    def __init__(self,particles,neighbour_list,cutoff=5.0):
        Force.__init__(self,particles,neighbour_list,cutoff=cutoff)

    def apply(self):
        """ Calculates spam interaction between all particles.
        """
        p = self.p
        nl = self.nl

        # Size parameters
        n = p.n
        d = p.dim
        ni = nl.nip
       
        # Particle properties
        x = np.reshape(p.r[0:n,:].copy(),(n,d),order='F')
        v = np.reshape(p.v[0:n,:].copy(),(n,d),order='F')
        a = np.zeros([n,d],order='F')
        T = p.t[0:n].reshape((n,1),order='F')
        rho = np.zeros((n),order='F')
        rho_lr = np.zeros((n),order='F')
        u = np.reshape(p.u[0:n].copy(),(n,1),order='F')
        sml = np.reshape(p.h[0:n].copy(),(n,1),order='F')
        sml_lr = np.reshape(p.hlr[0:n].copy(),(n,1),order='F')
        #jq = np.reshape(p.jq[0:n,:].copy(),(n,d)),order='F')
        jq = np.zeros([n,d],order='F')
        grad_rho = np.zeros((n,d),order='F')
        grad_v = np.zeros((n,d,d),order='F')
        grad_rho_lr = np.zeros((n,d),order='F')
        mass = np.reshape(p.m[0:n].copy(),(n,1),order='F')

        # Neighbourly properties
        ilist = nl.iap[0:ni,:].copy()+1
        dv =  np.reshape(nl.dv[0:ni,:].copy(),(ni,d),order='F')
        rij = np.reshape(nl.rij[0:ni].copy(),(ni,1),order='F')
        drij = np.reshape(nl.drij[0:ni,:].copy(),(ni,d),order='F')
        #w = np.zeros((ni),order='F') 
        #dwdx = np.zeros((ni,d),order='F') 
        #w_lr = np.zeros((ni),order='F') 
        #dwdx_lr = np.zeros((ni,d),order='F') 
      
        # Velocity diff
        sphlib.sphlib.calc_dv(dv,ilist,v)

        # Kernels and kernel gradients
        w,dwdx = fkernel.kernel.smoothing_kernels(rij,drij,ilist,sml,2)
        w_lr,dwdx_lr = fkernel.kernel.smoothing_kernels(rij,drij,ilist,sml_lr,2)

        # Density summation
        fkernel.kernel.density_sum(rho,grad_rho,ilist,sml,mass,w,dwdx,2)
        fkernel.kernel.density_sum(rho_lr,grad_rho_lr,ilist,sml_lr,mass,w_lr,dwdx_lr,2)

        # Pressure tensor components
        p_rev = np.zeros((n,d,d),order='F')
        p_rev_lr = np.zeros((n,d,d),order='F')
        pi_irr = np.zeros((n,d,d),order='F')

        # Speed of sound
        c = np.ones([n],dtype=float,order='F')
        
        # Constants
        eta = np.zeros([n],order='F')
        zeta = np.zeros([n],order='F')
        eta[:] = 1.0
        zeta[:] = 0.1

        dedt = np.zeros([n],dtype=float,order='F')
        
        feos.eos.calc_vdw_temp(u,T,rho)
        #feos.eos.calc_vdw_energy(u,T,rho)
        T [T < 0.0] = 0.0
        # print a.flags.f_contiguous

        # Call the force subroutine
        t1 = time()
        sphforce3d.sphforce3d.calc_sphforce3d( 
            ilist,x,v,a,  
            p_rev,p_rev_lr,pi_irr,          
            grad_rho_lr,grad_v,               
            u,dedt,mass,rho,rho_lr,T,jq,     
            c,eta,zeta,               
            dv,rij,drij,  
            sml,sml_lr,w,dwdx,dwdx_lr)#,n,ni)

        #np.set_printoptions(precision=5,suppress=True)
        feos.eos.calc_vdw_temp(u,T,rho)

        # Resend data to python object
        p.rho[0:n] = rho[0:n]
        p.rho_lr[0:n] = rho_lr[0:n]
        p.udot[0:n] = dedt[0:n]
        nl.wij[0:ni] = w[0:ni]
        nl.dwij[0:ni,:] = dwdx[0:ni,:]
        nl.wij_lr[0:ni] = w_lr[0:ni]
        nl.dwij_lr[0:ni,:] = dwdx_lr[0:ni,:]
        p.P[0:n] = p_rev + p_rev_lr + pi_irr
        #p.p = 
        #p.pco = 
        p.vdot[0:n,:] = a[0:n,:]
        #p.pco[0:n] = pco[0:n]
        p.t[0:n] = T[0:n,0]
        p.jq[0:n,:] = jq[0:n,:]

