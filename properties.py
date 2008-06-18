""" The idea behind this module - it is for calculating sph densities,
    interparticle distances and any other properties that are
    needed before computing the rates of change.
"""

import sys
sys.path.append("../vasp")
import spdensity
import spkernel
import math

H = 5.0
ADASH = 2.0
BDASH = 0.5
KBDASH = 1.0
RHONAUGHT = 1.0

def ideal_isothermal(rho):
    """ Calculates the pressure from the kinetic
        equation of state. 
        Isothermal equation of state
    """
    return (rho * KBDASH)


def art_water(rho):
    """ Equation of state. Isothermal, with a reference density.
        Models a compressible liquid.
    """
    return ((rho - RHONAUGHT)*KBDASH)


#calc_pressure = art_water
calc_pressure = ideal_isothermal

def spam_properties(p,nl):
    """ Calculates and assigns:
        interparticle distances and displacements
        kernel values
        kernel gradient values
        and smoothed particle
        summation densities for the particle data
        structure
    
        todo: move the spam stuff to sp_neighbour list
    
    """
    # self contribution to rho
    zerokern = spkernel.lucy_kernel(0.0,(0.0,0.0),H)[0]
    #for i in range(p.n):
    p.rho[0:p.n] = zerokern

    # calc the distances, kernels, densities
    for k in range(nl.nip):
        i = nl.iap[k,0]
        j = nl.iap[k,1]
        
        nl.drij[k,0] = p.r[j,0] - p.r[i,0]
        nl.drij[k,1] = p.r[j,1] - p.r[i,1]
        rsquared = nl.drij[k,0]**2 + nl.drij[k,1]**2 
        nl.rij[k] = math.sqrt(rsquared)

        nl.wij[k], nl.dwij[k] = spkernel.lucy_kernel(nl.rij[k],nl.drij[k,:],H)

        p.rho[i] += nl.wij[k] * p.m[j]
        p.rho[j] += nl.wij[k] * p.m[i]

    for i in range(p.n):
        p.p[i] = calc_pressure(p.rho[i])

        




