#! /usr/local/bin/python

""" 
    Batch script for Smooth Particle solver. 
    Copyright Andrew Charles 2009
    All rights reserved.
"""

import sys
from time import time
import particles
#import c_forces as forces
import forces
import neighbour_list
from properties import spam_properties
from spam_nc import create_sph_ncfile, write_step
import spam_complete_force

# Global variables
MAX_STEPS = 1000
NDIM = 3
XMAX = 7 
YMAX = 7
ZMAX = 7
VMAX = 0.0
dt = 0.05
SPACING = 0.9
LIVE_VIEW = False
SIDE = (5,5,5)
NP = SIDE[0] * SIDE[1] * SIDE[2]
TEMPERATURE = 1.0
HLONG = 4.0
HSHORT = 2.0

ofname = 'output.nc'

p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit='grid',vmax=VMAX
    ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX)

def initialise():
    global p
    print "Initialising"
    p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit='grid',vmax=VMAX
        ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX
        ,temperature=TEMPERATURE,hlong=HLONG,hshort=HSHORT)
    nl = neighbour_list.SortedVerletList(p,cutoff=4.0)
    p.nlists.append(nl)
    p.nl_default = nl
    p.forces.append(spam_complete_force.SpamComplete(p,nl))
    p.forces.append(forces.FortranCollisionForce(p,nl))
    nl.build()
    nl.separations()
    spam_properties(p,nl)
    cnt = 0
    attribs = {'name':'Andrew', 'age':33}
    create_sph_ncfile(ofname,attribs,NP,NDIM)

if __name__ == "__main__":
    initialise()
    print "STEP   INT  DERIV =  PAIR + SPAM +  FORCE   "
    for i in range(MAX_STEPS):
        tstart = time()
        p.update(dt)
        write_step(ofname,p)
        print "%5.3f  " %(time() - tstart)             \
             + "%5.3f  " %(p.timing['integrate time']) \
             + "%5.3f  " %(p.timing['deriv time'])     \
             + "%5.3f  " %(p.timing['pairsep time'])   \
             + "%5.3f  " %(p.timing['SPAM time'])     \
             + "%5.3f  " %(p.timing['force time'])
    


