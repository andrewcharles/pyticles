#! /usr/local/bin/python

""" 
    Batch script for Smooth Particle solver. 
    Copyright Andrew Charles 2009
    All rights reserved.
"""

import sys
from time import time
import particles
import neighbour_list
from properties import spam_properties
import spam_complete_force

# Global variables
MAX_STEPS = 5
NDIM = 3
XMAX = 12
YMAX = 12
ZMAX = 12
VMAX = 0.0
dt = 0.05
SPACING = 0.9
LIVE_VIEW = False
SIDE = (5,5,5)
NP = SIDE[0] * SIDE[1] * SIDE[2]
TEMPERATURE = 1.8
HLONG = 4.0
HSHORT = 2.0

p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit='grid',vmax=VMAX
    ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX
    ,temperature=TEMPERATURE,hlong=HLONG,hshort=HSHORT)

nl = neighbour_list.SortedVerletList(p,cutoff=4.0)
p.nlists.append(nl)
p.nl_default = nl
p.forces.append(spam_complete_force.SpamComplete(p,nl))
nl.build()
nl.separations()
spam_properties(p,nl)

print "STEP   INT  DERIV =  PAIR + SPAM +  FORCE   "
for i in range(MAX_STEPS):
    tstart = time()
    p.update(dt)
    print "%5.3f  " %(time() - tstart)             \
         + "%5.3f  " %(p.timing['integrate time']) \
         + "%5.3f  " %(p.timing['deriv time'])     \
         + "%5.3f  " %(p.timing['pairsep time'])   \
         + "%5.3f  " %(p.timing['SPAM time'])     \
         + "%5.3f  " %(p.timing['force time'])


