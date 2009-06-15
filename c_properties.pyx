""" The idea behind this module - it is for calculating sph densities,
    interparticle distances and any other properties that are
    needed before computing the rates of change.
    
    Copyright Andrew Charles 2008
    All rights reserved.
    This module is new BSD licensed.

"""

import sys
#sys.path.append("../vasp")
import spdensity
import spkernel
import math
import numpy as np
cimport numpy as np

H = 3.0
ADASH = 2.0
BDASH = 0.5
KBDASH = 1.0
RHONAUGHT = 1.0
ADKE = False
pi = math.pi

def ideal_isothermal(rho,t):
    """ Calculates the pressure from the kinetic
        equation of state. 
        Isothermal equation of state
    """
    return (rho * KBDASH)


def art_water(rho,t):
    """ Equation of state. Isothermal, with a reference density.
        Models a compressible liquid.
    """
    return ((rho - RHONAUGHT)*KBDASH)

def vdw(rho,t):
    """ Van der Waals repulsive pressure.
    """
    return (rho*KBDASH*t)/(1-rho*BDASH), - ADASH*rho*rho


#calc_pressure = art_water
calc_pressure = vdw


cdef void lucy_kernel_3d(float r,float dx[3],float h,float* w, float dwdx[3]):
    """High performance sph kernel. Thanks Klaus Dolag from whose
       notes I got the normalisation for 3D. Yes, I have no time for
       algebra. I really should get on top of Sage.
    """
    
    cdef float q

    q = 105. / (pi * 16 * (h**3) )

    r = abs(r)
    if( r > h ):
        w[0] = 0
        dwdx[0] = 0
        dwdx[1] = 0
        dwdx[2] = 0
        return

    w[0] = q * (1 + 3.*r/h)*((1.-r/h))**3

    if(r == 0):
        dwdx[0] = 0.0
        dwdx[1] = 0.0
        dwdx[2] = 0.0
    else:
        dwdx[0] =  q * ( (-12./(h**4))*(r**3) + \
                    (24./(h**3))*(r**2) - (12.*r/(h**2)) )* dx[0]/r
        dwdx[1] =  q * ( (-12./(h**4))*(r**3) + \
                    (24./(h**3))*(r**2) - (12.*r/(h**2)) )* dx[1]/r
        dwdx[2] =  q * ( (-12./(h**4))*(r**3) + \
                    (24./(h**3))*(r**2) - (12.*r/(h**2)) )* dx[2]/r


def spam_properties(p,nl,h):
    """ Calculates and assigns:
        kernel values
        kernel gradient values
        and smoothed particle
        summation densities for the particle data
        structure
        todo: move the spam stuff to sp_neighbour list
    
    """
    n = p.n

    # self contribution to rho
    cdef float zerokern
    cdef float wk
    cdef float dwk[3]
    cdef float rk, rksq
    cdef float drk[3]
    cdef float dvk[3]
    cdef int i,j

    cdef np.ndarray[np.int_t,ndim=2,mode='c'] _iap
    cdef np.ndarray[np.float_t,ndim=2,mode='c'] _drij
    cdef np.ndarray[np.float_t,ndim=2,mode='c'] _r
    cdef np.ndarray[np.float_t,ndim=1,mode='c'] _rho
    cdef np.ndarray[np.float_t,ndim=2,mode='c'] _v
    cdef np.ndarray[np.float_t,ndim=3,mode='c'] _gradv
    cdef np.ndarray[np.float_t,ndim=1,mode='c'] _rij
    cdef np.ndarray[np.float_t,ndim=1,mode='c'] _wij
    cdef np.ndarray[np.float_t,ndim=2,mode='c'] _dwij
    cdef np.ndarray[np.float_t,ndim=1,mode='c'] _m

    _iap = nl.iap.astype(np.int)
    _drij = nl.drij.astype(np.float)
    _rij = nl.rij.astype(np.float)
    _wij = nl.wij.astype(np.float)
    _dwij = nl.dwij.astype(np.float)
    _r = p.r.astype(np.float)
    _v = p.v.astype(np.float)
    _gradv = p.gradv.astype(np.float)
    _rho = p.rho.astype(np.float)
    _vdot = p.vdot.astype(np.float) 
    _m = p.m.astype(np.float) 

    drk[0] = 0.0
    drk[1] = 0.0
    drk[2] = 0.0
    lucy_kernel_3d(0.0,drk,h,&zerokern,dwk)

    for i in range(p.n):
        _rho[i] = zerokern
        _gradv[i] = 0.0

    # calc the distances kernels and densities
    # we redo the distances here because the particles may
    # have moved but we didn't rebuild the list
    # this functin should be renamed neibourly properties
    # or something like that
    for k in range(nl.nip):
        i = _iap[k,0]
        j = _iap[k,1]

        for q in range(3):
            dvk[q] = _v[j,q] - _v[i,q]

        rksq = 0
        for q in range(3):
            drk[q] += _r[j,q] - _r[i,q]
            _drij[k,q] = drk[q]
            rksq += drk[q]**2
        rk = math.sqrt(rksq)

        _rij[k] = rk
        lucy_kernel_3d(rk,drk,h,&wk,dwk)
        
        _wij[k] = wk
        
        for q in range(3):
            _dwij[k,q] = dwk[q]

        _rho[i] += _wij[k] * _m[j]
        _rho[j] += _wij[k] * _m[i]

        # BIG ERROR ALERT
        # LOOK UP THE 3d EXPRESSION FOR GRADV!!!!
        # 

        _gradv[i,0] += (_m[j]/_rho[j])*dvk[0]*_dwij[k,0]
        _gradv[i,1] += (_m[j]/_rho[j])*dvk[1]*_dwij[k,1]
        _gradv[i,2] += (_m[j]/_rho[j])*dvk[2]*_dwij[k,2]
        
        _gradv[j,0] -= (_m[i]/_rho[i])*dvk[0]*_dwij[k,0]
        _gradv[j,1] -= (_m[i]/_rho[i])*dvk[1]*_dwij[k,1]
        _gradv[j,2] -= (_m[i]/_rho[i])*dvk[2]*_dwij[k,2]

    if ADKE:
        # We are using adaptive density kernel estimation
        # the density calculated above was just a pilot
        # the smoothing length above is the reference length
        KSC = 1.0
        SENS = 0.5
        rhoav = np.mean(p.rho)
        p.h = H * KSC * ((p.rho/rhoav)**SENS)
        
    for i in range(p.n):
        # todo add some logic to determine whether we have a one or two part
        # pressure
        p.p[i],p.pco[i] = calc_pressure(p.rho[i],p.t[i])

    # Now pass the arrays back
    nl.wij[:] = _wij[:]
    nl.rij[:] = _rij[:]
    nl.dwij[:,:] = _dwij[:,:]
    nl.drij[:,:] = _drij[:,:]
    p.rho[:] = _rho[:]
    p.gradv[:,:] = _gradv[:,:]




