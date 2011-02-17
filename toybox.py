#! /usr/local/bin/python

""" 
    A very small system used for functional testing.

    A 10x10x10 scaling length box
    Scaling length 2.81nm
    Scaling time 1ns
    Scaling mass 7520 au
    a = 7.45e+04
    b = 5.84e-01
    kb = 3.29e+04
    Number of particles = 6x6x6
    Target density = M / V
    V = 1000


    Particle mass 0.1386
    Temperature 0.95

    ## Copyright Andrew Charles, RMIT 2010                ###
    ## All rights reserved.                               ###

"""

import sys
from time import time
import particles
import spam_complete_force
import neighbour_list
from properties import spam_properties
from spam_nc import create_sph_ncfile, write_step
import numpy as np
import spkernel, spdensity


def spawn_system(
        max_steps = 25,
        xmax = 10,
        YMAX = 10,
        ZMAX = 10 ,
        NDIM = 3,
        SIDE = (5,5,5),
        VMAX = 0.0,
        dt = 0.05,
        SPACING = 1.0,
        TEMPERATURE = 0.95,
        HLONG = 4.0,
        HSHORT = 2.0,
        RINIT = 'grid',
        ascl = 7.45e+04,
        bscl = 5.84e-01,
        kbscl = 3.29e+04,
        pmass = 1.386e-01,
        ofname = 'data/toybox.nc'
        ):
        """ Create a smooth particle system and integrate it forward in time.

            Returns the path of the netcdf file containing the results.
        """

        NP = SIDE[0]*SIDE[1]*SIDE[2]
        cnt = 0
        fps = 0

        print "Initialising"
        p = particles.SmoothParticleSystem(
                NP,maxn=NP,
                d=3,
                rinit=RINIT,
                vmax=VMAX,
                side=SIDE,
                spacing=SPACING,
                xmax=xmax,
                ymax=YMAX,
                zmax=ZMAX,
                temperature=TEMPERATURE,
                hlong=HLONG,
                hshort=HSHORT,
                thermostat_temp=TEMPERATURE,
                thermostat=True,
                mass=pmass
            )
        nl = neighbour_list.VerletList(p,cutoff=HLONG)
        p.nlists.append(nl)
        p.nl_default = nl
        p.forces.append(
            spam_complete_force.SpamComplete(
                p,nl,adash=ascl,bdash=bscl,kbdash=kbscl))
        #p.forces.append(forces.FortranCollisionForce(p,nl,cutoff=0.5))
        tstart = time()
        nl.build()
        nl.separations()
        spam_properties(p,nl)
        print 'Built list and calc properties',time()-tstart
        cnt = 0
        attribs = {'creator':'Andrew', 'log':'functional test'}
        create_sph_ncfile(ofname,attribs,NP,NDIM)
        print "STEP   INT  DERIV =  PAIR + SPAM +  FORCE   "
        tstartrun = time()
        for i in range(max_steps):
            tstart = time()
            p.update(dt)
            if np.isnan(p.r).any():
                print 'stopping due to nan'
                break
            if i % 10 == 0:
                write_step(ofname,p)
            print 'Step',i,'took',time()-tstart
            g = p.timing.keys()
            g.sort()
            for k in g:
                print k,p.timing[k]
        print 'Completed',i,'steps, in',time()-tstartrun
        return ofname

if __name__ == '__main__':
    spawn_system()


