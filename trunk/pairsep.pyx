""" 
    The most expensive part of a neighbour list is the computation of pair
    separations.

    Andrew Charles
"""

import numpy as np
cimport numpy as np
import math
import scipy
import neighbour_list
import particles

# CONSTANTS

ctypedef np.float_t DTYPE_t
DTYPE = np.float

cdef extern from "math.h":
    float sinf(float theta)
    float sqrtf(float a)
    float powf(float a, float b)

def pairsep(nl):
    """ Computes the seperations between a set of points.
        r[n,3] = array of positions
        iap[nk,2] = integer array of interacting pairs
        dr[nk,3] = returned array of pair separations
        ds[nk] = returned array of pair separation distances

    """
        
    cdef np.ndarray[np.int_t,ndim=2,mode='c'] _pairs 
    cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _r
    cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _v
    cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _dr 
    cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _dv
    cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _ds
    cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _rsq
    cdef float drx, dry, drz
    cdef unsigned int i,j,k

    n = nl.particle.n
    nk = nl.nip
    dim = nl.particle.dim

    _pairs = nl.iap.astype(np.int)
    _r = nl.particle.r.astype(DTYPE)
    _v = nl.particle.v.astype(DTYPE)
    
    _dr = np.zeros((nk,dim))
    _ds = np.zeros(nk)
    _dv =  np.zeros((nk,dim))
    _rsq =  np.zeros(nk)

    # Which is faster?
    #_dr = nl.rij.astype(np.float)
    #_ds = nl.rij.astype(np.float)
    #_dv = nl.dv.astype(np.float)
    #_rsq = nl.rsq.astype(np.float)

    for k in xrange(nk):
        i = _pairs[k,0]
        j = _pairs[k,1]
        drx = _r[j,0] - _r[i,0]
        dry = _r[j,1] - _r[i,1]
        drz = _r[j,2] - _r[i,2]
        _rsq[k] = powf(drx,2) + powf(dry,2) + powf(drz,2)
        _ds[k] = sqrtf(_rsq[k])
        _dr[k,0] = drx
        _dr[k,1] = dry
        _dr[k,2] = drz
        _dv[k,0] = _v[j,0] - _v[i,0]
        _dv[k,1] = _v[j,1] - _v[i,1]
        _dv[k,2] = _v[j,2] - _v[i,2]

    nl.drij[0:nk,:] = _dr[0:nk,:]
    nl.rij[0:nk] = _ds[0:nk]
    nl.rsq[0:nk] = _rsq[0:nk]
    nl.dv[0:nk,:] = _dv[0:nk,:]

