#! /usr/local/bin/python

""" 
    Optimised smooth particle implementation using Fortran
    
    Attempt to use the parameters of Bedaux et al

    Scaling parameters:
    t_s = 562K
    rho_s = 3.12 x 10^4 molm^-3
    x_s = 2.81 nm
    p_s = 200 bar

    At 200 bar, 561.6K
    rho_molar_liq = 3.1 x 10^4 molm^-3
    rho_molar_gas = 0.6
    
    Molar mass of water is 18.0153 gmol-1

    Box size is

    36 x 12 x 12

    The initial plan was to have one 6 x 12 x 12 end 
    constrained to liquid and one 6 x 12 x 12 end 
    constrained to vapour. I am no longer doing this.

    Copyright Andrew Charles 2009
    All rights reserved.
"""

import sys
from time import time
import particles
import spam_complete_force
import profile
import neighbour_list
from pyglet.gl import *
from properties import spam_properties
import numpy as np
import spkernel, spdensity

# Global variables
MAX_STEPS = 1000

XMAX = 36 
YMAX = 12
ZMAX = 12 
SIDE = (18,6,6)
VMAX = 0.0
dt = 0.05
SPACING = 1.0
NP = SIDE[0]*SIDE[1]*SIDE[2]
TEMPERATURE = 1.0
HLONG = 4.0
HSHORT = 2.0
RINIT = 'grid'

p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit=RINIT,vmax=VMAX
    ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX)
s = pview.ZPRView(p)
nl = neighbour_list.VerletList(p,cutoff=4.0)

cnt = 0
fps = 0

ofname = 'output.nc'
print "Initialising"
p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit=RINIT,vmax=VMAX
,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX
,temperature=TEMPERATURE,hlong=HLONG,hshort=HSHORT,
thermostat_temp=TEMPERATURE,thermostat=True)
nl = neighbour_list.VerletList(p,cutoff=5.0)
p.nlists.append(nl)
p.nl_default = nl
p.forces.append(spam_complete_force.SpamComplete(p,nl))
#p.forces.append(forces.FortranCollisionForce(p,nl,cutoff=0.5))
nl.build()
nl.separations()
spam_properties(p,nl)
cnt = 0
attribs = {'name':'Andrew', 'age':33}
create_sph_ncfile(ofname,attribs,NP,NDIM)
initialise()
print "STEP   INT  DERIV =  PAIR + SPAM +  FORCE   "
for i in range(MAX_STEPS):
    tstart = time()
    p.update(dt)
    if np.isnan(p.r).any():
        print 'stopping due to nan'
        break
    if i % 10 == 0:
        write_step(ofname,p)
print 'Completed',i,'steps'






