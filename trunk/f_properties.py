""" Wrapper around Fortran SPH code to compute SPH properties """
import sys
sys.path.append('/Users/acharles/masters/active/fsph')
import fkernel
import sphlib
import feos
from eos import vdw_co, vdw_hc
import numpy as np
import properties

feos.eos.adash = 2.0
feos.eos.bdash = 0.5
feos.eos.kbdash = 1.0

def spam_properties(p,nls,nll,hs,hl):
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
        (although not long range estimates of density, just longer h in 
        the pco calculation).

        Still no heat flux

    """

    n = p.n
    d = p.dim
    ni = nls.nip
    nil = nll.nip

    v = np.reshape(p.v[0:n,:].copy(),(n,d))#,order='F')
    dv =  np.reshape(nls.dv[0:ni,:].copy(),(ni,d))#,order='F')
    #dv = reshape(nl.dv[0:ni,:].copy(),(ni,d),order='F')#.transpose()
    nlist_short = np.reshape(nls.iap[0:ni,:].copy(),(ni,2))+1
    nlist_long = np.reshape(nll.iap[0:nil,:].copy(),(nil,2))+1
    
    #nlist = np.reshape(nlist,(ni,2),order='F')
    #nlist = nlist.transpose()
    #dv = dv.transpose()
    #v = v.transpose()
    #nlist = np.array(nlist[:,:].transpose(),order='F')

    T = p.t[0:n].reshape((n))#,order='F')
    sml = np.reshape(p.h[0:n].copy(),(n,1))#,order='F')
    sml_lr = np.reshape(p.hlr[0:n].copy(),(n,1))#,order='F')

    # This will not be necessary if the lists are sorted
    # (having two sets of distances)
    rij = np.reshape(nls.rij[0:ni].copy(),(ni,1))#,order='F')
    drij = np.reshape(nls.drij[0:ni,:].copy(),(ni,d))#,order='F')
   
    rij_lr = np.reshape(nll.rij[0:nil].copy(),(nil,1))#,order='F')
    drij_lr = np.reshape(nll.drij[0:nil,:].copy(),(nil,d))#,order='F')

    w = np.zeros((ni),order='F') 
    dwdx = np.zeros((ni,d),order='F') 
   
    w_lr = np.zeros((nil),order='F') 
    dwdx_lr = np.zeros((nil,d),order='F') 

    rho = np.zeros((n),order='F')
    u = np.reshape(p.u[0:n].copy(),(n,1))#,order='F')
    gradv = np.reshape(p.gradv[0:n,:,:].copy(),(n,d,d))#,order='F')
    grad_rho = np.zeros((n,d))#,order='F')
    mass = np.reshape(p.m[0:n].copy(),(n,1))#,order='F')

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
    sphlib.sphlib.calc_dv(dv,nlist_short,v)

    # Kernels and kernel gradients
    w,dwdx = fkernel.kernel.smoothing_kernels(rij[0:ni],drij[0:ni,:] \
        ,nlist_short,sml,1)
    w_lr,dwdx_lr = fkernel.kernel.smoothing_kernels(rij_lr,drij_lr \
        ,nlist_long,sml_lr,1)

    # Density summation
    fkernel.kernel.density_sum(rho,grad_rho,nlist_short,sml,mass,w,dwdx,1)
    print rho[0]

    # Grad v
    grad_v = np.zeros((n,d,d),order='F')
    sphlib.sphlib.calc_grad_v(grad_v,nlist_short,dwdx,dv,mass,rho)

    # (Python implementation)
    phc = p.p[0:n]
    pco = p.pco[0:n]
  
    print 'b4'
    print 'frho',rho[0]
    print u[0]
    print T[0]
    feos.eos.calc_vdw_temp(u,T,rho)
    print 'afta2'
    print u[0]
    print T[0]
    print properties.vdw_temp(rho[0],u[0])
    print properties.vdw_energy(rho[0],1.2)
    1/0


    feos.eos.calc_vdw_energy(u,T,rho)
    feos.eos.calc_vdw_temp(u,T,rho)
    feos.eos.calc_vdw_energy(u,T,rho)

    feos.eos.calc_vdw_hc_pressure(phc,rho,u)
    feos.eos.calc_vdw_cohesive_pressure(pco,rho,u)

    #phc = vdw_hc(rho,T)
    #pco = vdw_co(rho,T)

    # Resend data to python object
    p.rho[0:n] = rho[0:n]
    nls.wij[0:ni] = w[0:ni]
    nls.dwij[0:ni,:] = dwdx[0:ni,:]
    nll.wij[0:nil] = w_lr[0:nil]
    nll.dwij[0:nil,:] = dwdx_lr[0:nil,:]
    p.p[0:n] = phc[0:n]
    p.pco[0:n] = pco[0:n]
    p.t[0:n] = T[0:n]
    






