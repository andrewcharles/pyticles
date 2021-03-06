""" Wrapper around Fortran SPH code to compute SPH properties """
import sys
sys.path.append('/Users/acharles/masters/active/fsph')
import fkernel
import splib
import feos
#from eos import vdw_co, vdw_hc
import numpy as np
import properties

feos.eos.adash = 2.0
feos.eos.bdash = 0.5
feos.eos.kbdash = 1.0
thermalk = 1.0

def spam_properties(p,nl,hs,hl):
    """ Calculates and assigns:

        * kernel values
        * kernel gradient values
        * smoothed particle summation densities
        * velocity gradient

        Particle arrays will not be in fortran order.
        See if we can just do copies ($$$$$).
        Despite the copies it's still loads faster than straight python
        and seems to have an advantage over cython.

        nl.iap -- integer list of interactions
        nl.dv -- velocity difference

        p.rho -- particle density
        p.gradv -- particle velocity gradient
       
        I'm writing this to work with the long and short smoothing length

        Doing:
        * Heat flux

        To Do:
        * Long range density
        * Full short range reversible and irreversible tensors
        * Full long range reversible and irreversible tensors

    """

    n = p.n
    d = p.dim
    ni = nl.nip

    # If an array is input/output, it should have the array order set to fortran
    # v = np.reshape(p.v[0:n,:].copy(),(n,d))#,order='F')

    v = np.reshape(p.v[0:n,:].copy(),(n,d))
    dv =  np.reshape(nl.dv[0:ni,:].copy(),(ni,d))
    nlist = nl.iap[0:ni,:].copy()+1
    
    #nlist = nlist.transpose()
    #dv = dv.transpose()
    #v = v.transpose()
    #nlist = np.array(nlist[:,:].transpose(),order='F')

    T = p.t[0:n].reshape((n))#,order='F')
    sml = np.reshape(p.h[0:n].copy(),(n,1))#,order='F')
    sml_lr = np.reshape(p.hlr[0:n].copy(),(n,1))#,order='F')

    rij = np.reshape(nl.rij[0:ni].copy(),(ni,1))#,order='F')
    drij = np.reshape(nl.drij[0:ni,:].copy(),(ni,d))#,order='F')
   
    w = np.zeros((ni),order='F') 
    dwdx = np.zeros((ni,d),order='F') 
   
    w_lr = np.zeros((ni),order='F') 
    dwdx_lr = np.zeros((ni,d),order='F') 

    rho = np.zeros((n))#,order='F')
    rho_lr = np.zeros((n))#,order='F')
    u = np.reshape(p.u[0:n].copy(),(n,1))#,order='F')
    gradv = np.reshape(p.gradv[0:n,:,:].copy(),(n,d,d))#,order='F')
    jq = np.reshape(p.jq[0:n,:].copy(),(n,d))#,order='F')
    grad_rho = np.zeros((n,d))
    grad_rho_lr = np.zeros((n,d))
    mass = np.reshape(p.m[0:n].copy(),(n,1))

    if ni == 0:
        # set values to isolated particle values
        # and issue a warning
        print 'Warning! No interactions!'
        return

    #print nlist.flags

    # Velocity diff
    #dv = sphlib.sphlib.calc_dv(dv,nlist,v)
    #print dv.flags
    #print v.flags
    splib.splib.calc_dv(dv,nlist,v)

    # Kernels and kernel gradients
    w,dwdx = fkernel.kernel.smoothing_kernels(rij[0:ni],drij[0:ni,:] \
        ,nlist,sml,2)
    w_lr,dwdx_lr = fkernel.kernel.smoothing_kernels(rij[0:ni],drij[0:ni] \
        ,nlist,sml_lr,2)
   
    # Density summation
    fkernel.kernel.density_sum(rho,grad_rho,nlist,sml,mass,w,dwdx,2)
    fkernel.kernel.density_sum(rho_lr,grad_rho_lr,nlist,sml_lr,mass,w_lr,dwdx_lr,2)

    # Grad v
    grad_v = np.zeros((n,d,d),order='F')
    splib.splib.calc_grad_v(grad_v,nlist,dwdx,dv,mass,rho)

    phc = np.reshape(p.p[0:n].copy(),(n),order='F')
    pco = np.reshape(p.pco[0:n].copy(),(n),order='F')
  
    feos.eos.calc_vdw_energy(u,T,rho)
    T [T < 0.0] = 0.0
    feos.eos.calc_vdw_hc_pressure(phc,rho,u)
    feos.eos.calc_vdw_cohesive_pressure(pco,rho_lr,u)
    feos.eos.calc_vdw_temp(u,T,rho)

    # Capillary pressure
     
    # Heat Flux
    # No idea why this needs to be order F, and can't be the copy of the
    # particle's jq.
    jq = np.zeros((n,d),order='F')
    
    # subroutine calc_heat_flux_3d(q,ilist,rho,m,tmp,dwdx_jij,thermalk,n,ni)
    splib.splib.calc_heat_flux_3d(jq,nlist,rho,mass,T,-dwdx,thermalk)

    # Python implementation
    #phc = vdw_hc(rho,T)
    #pco = vdw_co(rho,T)
    # viscosity

    # Resend data to python object
    p.rho[0:n] = rho[0:n]
    p.rho_lr[0:n] = rho_lr[0:n]
    nl.wij[0:ni] = w[0:ni]
    nl.dwij[0:ni,:] = dwdx[0:ni,:]
    nl.wij_lr[0:ni] = w_lr[0:ni]
    nl.dwij_lr[0:ni,:] = dwdx_lr[0:ni,:]
    p.p[0:n] = phc[0:n]
    p.pco[0:n] = pco[0:n]
    p.t[0:n] = T[0:n]
    p.jq[0:n,:] = jq[0:n,:]
    






