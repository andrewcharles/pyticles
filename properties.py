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
import numpy

H = 3.0
ADASH = 2.0
BDASH = 0.5
KBDASH = 1.0
RHONAUGHT = 1.0
ADKE = False

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

def spam_properties(p,nl,h):
    """ Calculates and assigns:

        * kernel values
        * kernel gradient values
        * smoothed particle summation densities
        * velocity gradient
    
    """
    # self contribution to density
    zerokern = spkernel.lucy_kernel(0.0,(0.0,0.0),h)[0]
    p.rho[0:p.n] = zerokern
    p.gradv[0:p.n] = 0.0

    # calc the kernels, velocities and densities
    for k in range(nl.nip):
        i = nl.iap[k,0]
        j = nl.iap[k,1]

        nl.wij[k], nl.dwij[k] = spkernel.lucy_kernel(nl.rij[k],nl.drij[k,:],h)

        p.rho[i] += nl.wij[k] * p.m[j]
        p.rho[j] += nl.wij[k] * p.m[i]
    
        dv = p.v[j,:] - p.v[i,:]

        for a in range(p.dim):
            for b in range(p.dim):
                p.gradv[i,a,b] += (p.m[j]/p.rho[i])*dv[a]*nl.dwij[k,b]
                p.gradv[j,a,b] += (p.m[i]/p.rho[j])*dv[a]*nl.dwij[k,b]

    if ADKE:
        # We are using adaptive density kernel estimation
        # the density calculated above was just a pilot
        # the smoothing length above is the reference length
        KSC = 1.0
        SENS = 0.5
        rhoav = numpy.mean(p.rho)
        p.h = H * KSC * ((p.rho/rhoav)**SENS)
        

    for i in range(p.n):
        # todo add some logic to determine whether we have a one or two part
        # pressure
        p.p[i],p.pco[i] = calc_pressure(p.rho[i],p.t[i])

        




