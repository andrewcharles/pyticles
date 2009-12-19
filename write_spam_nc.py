#! /usr/local/bin/python

""" 
    Create and write a netcdf file of SPAM output.

    todo - would be nice to have a dimension for every timestep, and a dimension
    just for the frames we've snapshotted. should read up on how to do this beofre
    making changes that might cause regression.
"""

import numpy
import netCDF4 
import sphstate
import os
import glob
import sys
import load_sphvars


def create_sph_ncfile(filename,attribs,n,dim):
    """ Create a netCDF file for SPAM output data, roughly
        following the AMBER specification.
        see http://amber.scripps.edu/netcdf/nctraj.html  for inspiration
        and some of the constants(?)

        attribs -- a dictionary of attributes for the nc file name, value
        n -- number of particles
        dim -- spatial dimensions

    """
    nc_file = netCDF4.Dataset(filename,'w')

    # Miscellaneous attributes
    setattr(nc_file,'Date',1)
    setattr(nc_file,'Creator','ac')

    # sphvars file attributes
    for name,val in attribs.iteritems():
        setattr(nc_file,name,val)
 
    # Create netcdf dimensions
    # number of particles
    # spatial dimensions
    # timestep number
    nc_file.createDimension('timestep',None)
    nc_file.createDimension('particle',n)
    nc_file.createDimension('spatial',dim)

    # Create variables for the dimensions, and populate them
    tstep = nc_file.createVariable('timestep','d',('timestep',))
    part = nc_file.createVariable('particle','i',('particle',))
    space = nc_file.createVariable('spatial','i',('spatial',))
   
    part[:] = numpy.array(range(n))
    space[:] = numpy.array([0,1,2])

    dimnames = nc_file.dimensions.keys()

    # Set up variables
    # every particle property has a variable
    # and there are also variables for the box size
    # and the box dimensions.
    # a variable for 'time elapsed' at each step (for variable stepping)

    #each variable needs a "units" attribute

    #vector variables
    v_dims =('timestep','particle','spatial')

    #scalar variables
    sc_dims = ('timestep','particle')
    
    #histogram variables
    hist_dims = ('timestep')

    #total and average variables
    tot_dims = ('timestep')

    r = nc_file.createVariable('position','d',v_dims)
    v = nc_file.createVariable('velocity','d',v_dims)
    #a = nc_file.createVariable('acceleration','d',v_dims)
    #temp = nc_file.createVariable('temperature','d',sc_dims)
    energy = nc_file.createVariable('internal_energy','d',sc_dims)
    mass = nc_file.createVariable('mass','d',sc_dims)
    #rho = nc_file.createVariable('density','d',sc_dims)
    #press = nc_file.createVariable('pressure','d',sc_dims)
    #ss =nc_file.createVariable('sound_speed','d',sc_dims)
    #visc =nc_file.createVariable('viscosity','d',sc_dims)
    #h = nc_file.createVariable('smoothing_length','d',sc_dims)
    #hl = nc_file.createVariable('long_smoothing_length','d',sc_dims)
    #q = nc_file.createVariable('heat_flux','d',v_dims)
    #vsm= nc_file.createVariable('smoothed_velocity','d',v_dims)
    #psm =nc_file.createVariable('smoothed_pressure','d',sc_dims)
    #tmpsm =nc_file.createVariable('smoothed_temperature','d',sc_dims)
    #grad_rho = nc_file.createVariable('density_gradient','d',v_dims)
    #ptype = nc_file.createVariable('particle_type','u1',sc_dims)

    #now set up the non-particle averaged or total system variables
    # kinetic energy, internal energy, isolated Hamiltonian

    #V = nc_file.createVariable('total_kinetic_energy','d',tot_dims) 
    #T = nc_file.createVariable('total_internal_energy','d',tot_dims)
    #tav = nc_file.createVariable('average_temp','d',tot_dims)
    #rhoav = nc_file.createVariable('rho_average','d',tot_dims)
    #tstat_energy = nc_file.createVariable('thermostat_energy','d',tot_dims)
    #TV =  nc_file.createVariable('hamiltonian','d',tot_dims)
    #dti = nc_file.createVariable('dt','d',tot_dims)
    #sys_dt = nc_file.createVariable('systime','d',tot_dims)
     
    nc_file.close()

def write_step(filename,p):
    """ Write out this step to the output file. 
    
        p -- particle object
    
    """
    
    ncfile = netCDF4.Dataset(filename,'r+')
    time = ncfile.dimensions['timestep']
    i = len(time)
    print i,'timesteps'
    timev = ncfile.variables['timestep']
    timev[i] = i+1

    r = ncfile.variables['position']
    v = ncfile.variables['velocity']
    u = ncfile.variables['internal_energy']
    m = ncfile.variables['mass']
    n = p.n
        
    r[i,:,:] = p.r[0:n,:]
    v[i,:,:] = p.v[0:n,:]
    u[i,:] = p.u[0:n]
    m[i,:] = p.m[0:n]

    ncfile.sync()
    ncfile.close()

if __name__ == "__main__":
    print 'Testing'
    import particles

    # TEST NC WRITE
    # 1. Create a netcdf file
    attribs = {'name':'Andrew', 'age':33}
    filename = 'test.nc'
    create_sph_ncfile(filename,attribs,27,3)

    # 2. Create a particle
    SIDE = (3,3,3)
    p = particles.SmoothParticleSystem(27,maxn=30,d=3,rinit='grid',vmax=1.0
        ,side=SIDE,spacing=1.0,xmax=10,ymax=10,zmax=10)

    # 3. Write the particle data out
    write_step(filename,p)

