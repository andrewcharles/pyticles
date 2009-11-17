#!/usr/local/bin/python
""" Renders a spam field in two dimensions.
    Copyright Andrew Charles 2008
    All rights reserved.
    This module is new BSD licensed.
"""

import math
import numpy
import spkernel
import interpolate
import neighbours
import pylab
import configuration

NPTS = 100
# resolution required
nx = 0
ny = 0
ngridpts = 0

def map_to_grid():
    """ Maps particle coordinates to grid coordinates """
    print "wtf"

def get_grid_map(xmin,xmax,ymin,ymax,res):
    """ Returns the coordinates of the centers of the
        grid points.
        res is the resolution of the grid, i.e. the spacing
    """
    #global nx
    #global ny
    #global ngridpts
    gridx = numpy.arange(xmin+res/2.0,xmax,res)
    gridy = numpy.arange(ymin+res/2.0,ymax,res)
    #gridx = numpy.arange(0+res/2.0,xmax,res)
    #gridy = numpy.arange(0+res/2.0,ymax,res)
    #nx = int(((xmax - xmin)/res) )
    #ny = int(((ymax - ymin)/res) )
    #ngridpts = nx*ny

    return gridx,gridy 

def func2(x,y):
    return numpy.exp(x+y)


def renderspam2d(x,y,r,m,rho,a,h,xmin,xmax,ymin,ymax,cutoff):
    """ x: x coordinates of the grid
        y: y coordinates of the grid
        r: coordinates of points [x,y]
        rho: sph summation densities of the points
        m: mass of each point
        a: property we are rendering (of each point) 
        xmin,xmax,ymin,ymax
    """

    bounds = xmin,xmax,ymin,ymax
   
#    Z = interpolate.map_to_grid(x,y,r,m,rho,a,h,bounds,cutoff) 
    Z = interpolate.splash_map(x,y,r,m,rho,a,h,bounds,cutoff) 

    # problem with imshow is that it does its
    # own interpolation - I want to use the sph
    # kernel only
    # No idead why the transpose is needed
    # It appears that contour and imshow reverse the usual row/column order
    Z = Z.transpose()
    #pylab.contour(x,y,Z)
    #pylab.colorbar()
    pylab.axis([xmin,xmax,ymin,ymax])
    pylab.imshow( Z , origin='lower', extent=(xmin,xmax,ymin,ymax),interpolation='bilinear',vmin=0.0,vmax=2.0 )
    #pylab.imshow( Z , origin='lower', extent=(xmin,xmax,ymin,ymax),interpolation='bilinear')
    pylab.colorbar()
    pylab.axis([xmin,xmax,ymin,ymax])


#Random plot stuff
#        axis([0,xmax,0,ymax])
#        rhoxy = array( rho[i,:])
#        contour(r[i,:,0],r[i,:,1],rhoxy)
    #pylab.imshow( Z , interpolation='bilinear', origin='lower' )
        #contourf(Z)


def main():
    """ Plots a regular 2d spam field
    """

    # initial positions
    xmin = 0
    xmax = 5
    ymin = 0
    ymax = 10
    n = 10
    cutoff = 3
    r = configuration.random(n,xmin,xmax,ymin,ymax)

    h = numpy.zeros((r.shape[0]))
    h[:] = 2
    #ry = numpy.array([1,1,1,2,2,2,3,3,3])
    gridx,gridy = get_grid_map(xmin,xmax,ymin,ymax)
    pylab.ion()
    pylab.axis([xmin,xmax,ymin,ymax])
    renderspam2d(gridx,gridy,r,h,xmin,xmax,ymin,ymax,cutoff)
    pylab.plot(r[:,0],r[:,1],'wo')
    pylab.axis([xmin,xmax,ymin,ymax])
    pylab.show()

if __name__ == "__main__":
    main()

