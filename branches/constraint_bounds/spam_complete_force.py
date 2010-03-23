from forces import Force

""" A wrapper around the fortran sphforce3d routine.
    Computes some smooth particle properties first.

    #np.set_printoptions(precision=5,suppress=True)

"""

import sys
sys.path.append('/Users/acharles/masters/active/fsph')
import fkernel
import splib
import feos
import numpy as np
import sphforce3d
from time import time


class SpamComplete(Force):
    """ Compute all the smooth particle properties and forces
        in a single function that uses the fortran routines.

        Parameters:
            adash -- van der Waals attraction (sets in feos module)
            bdash -- van der Waals repulsion (sets in feos module)
            kbdash -- van der Waals k (sets in feos module)
            sigma -- repulsive core size (fortran module variable)
            rcoef -- repulsive core strength (fortran module variable)
            cgrad -- density gradient coefficient
            eta -- shear viscosity
            zeta -- bulk viscosity
            kernel_type -- 1: Gauss, 2: Lucy, 3: Debrun

    """

    def __init__(self,particles,neighbour_list,
        adash=2.0,
        bdash=0.5,
        kbdash=1.0,
        sigma=0.0,
        rcoef=0.0,
        cgrad=1.0,
        eta=1.0,
        zeta=0.1,
        kernel_type=2,
        cutoff=5.0):

        # Call the generic Force init method
        Force.__init__(self,particles,neighbour_list,cutoff=cutoff)

        # Set the parameters in this module and imported modules
        feos.eos.adash = adash
        feos.eos.bdash = bdash
        feos.eos.kbdash = kbdash
        sphforce3d.sigma = sigma
        sphforce3d.rcoef = rcoef
        sphforce3d.cgrad = cgrad
        self.eta = eta
        self.zeta = zeta
        self.kernel_type = 2 


    def apply(self):
        """ Calculates spam interaction between all particles.
        """
        p = self.p
        nl = self.nl

        # Size parameters
        n = p.n
        d = p.dim
        ni = nl.nip
       
        # Copy / reference particle properties
        x = np.asfortranarray(p.r[0:n,:])
        v = np.asfortranarray(p.v[0:n,:])
        T = np.asfortranarray(p.t[0:n])
        mass = np.asfortranarray(p.m[0:n])
        rho = np.zeros((n),order='F')
        rho_lr = np.zeros((n),order='F')
        u = np.asfortranarray(p.u[0:n])
        sml = np.asfortranarray(p.h[0:n])
        sml_lr = np.asfortranarray(p.hlr[0:n])
        
        # Copy / reference Neighbourly properties
        ilist = np.asfortranarray((nl.iap[0:ni,:]+1))
        dv =  np.asfortranarray(nl.dv[0:ni,:])
        rij = np.asfortranarray(nl.rij[0:ni])
        drij = np.asfortranarray(nl.drij[0:ni,:])
        
        # New allocations
        a = np.zeros([n,d],order='F')
        #x = np.zeros([n,d],order='F')
        jq = np.zeros([n,d],order='F')
        grad_rho = np.zeros((n,d),order='F')
        grad_v = np.zeros((n,d,d),order='F')
        grad_rho_lr = np.zeros((n,d),order='F')

        # Velocity diff
        splib.splib.calc_dv(dv,ilist,v)

        # Kernels and kernel gradients
        w,dwdx = fkernel.kernel.smoothing_kernels(rij,drij,ilist,sml,
            self.kernel_type)
        w_lr,dwdx_lr = fkernel.kernel.smoothing_kernels(rij,drij,ilist,sml_lr,
            self.kernel_type)

        # Density summation
        fkernel.kernel.density_sum(rho,grad_rho,ilist,sml,mass,w,dwdx,
            self.kernel_type)
        fkernel.kernel.density_sum(rho_lr,grad_rho_lr,ilist,sml_lr,mass,w_lr,
            dwdx_lr,self.kernel_type)

        feos.eos.calc_vdw_temp(u,T,rho)

        # Pressure tensor components
        p_rev = np.zeros((n,d,d),order='F')
        p_rev_lr = np.zeros((n,d,d),order='F')
        pi_irr = np.zeros((n,d,d),order='F')

        # Speed of sound
        c = np.ones([n],dtype=float,order='F')
        
        # Constants
        eta = np.zeros([n],order='F')
        zeta = np.zeros([n],order='F')
        eta[:] = self.eta 
        zeta[:] = self.zeta

        dedt = np.zeros([n],dtype=float,order='F')
        feos.eos.calc_vdw_temp(u,T,rho)
        T [T < 0.0] = 0.0
        # print a.flags.f_contiguous

        # Call the force subroutine
        sphforce3d.sphforce3d.calc_sphforce3d( 
            ilist,x,v,a,  
            p_rev,p_rev_lr,pi_irr,          
            grad_rho_lr,grad_v,               
            u,dedt,mass,rho,rho_lr,T,jq,     
            c,eta,zeta,               
            dv,rij,drij,  
            sml,sml_lr,w,dwdx,dwdx_lr)
        
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
        p.p[0:n] = p_rev[:,0,0]
        p.pco[0:n] = p_rev_lr[:,0,0]
        p.vdot[0:n,:] = a[0:n,:]
        p.t[0:n] = T[0:n]
        p.jq[0:n,:] = jq[0:n,:]

