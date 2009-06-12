""" Part of the sph suite
    this module performs an sph interpolation
    over a disordered set of points onto a
    regular grid.

    Big up to Daniel Price and his code for giving me the idea to
    loop over particles to determine which pixels they contribute to
    before building the list.
    Price (2007), Publ. Astron. Soc. Aust., 24, 159-173 (arXiv:0709.0832) 
    
    As DPS's code is is GPLed, so is this module.

"""
from math import *
import numpy
import neighbours
import spkernel

def smooth_point(r,dr,x,vol,h,type):
    """ Calculates the smoothed value at a point given a set of
        distances and property values.
        r[n]: the distance between this point and all neighbours
        x[n]: the property value at each point
        vol[n]: mj / rhoj, volumes of particles
        h[n]: array of smoothing lengths
        type: string, e.g "lucy"
        <A> = sum ( (vj) * Aj )
    """

    # Should specify this as a double
    x_smooth = 0.0
    j = 0
    for rj in r:
        x_smooth = x_smooth + x[j] * vol[j] * spkernel.kernel(rj,dr[j],h[j],type)[0]
        j += 1
    return x_smooth

def splash_map(x,y,r,m,rho,a,h,bounds,cutoff): 
    """ I used Daniel Price's algorithm for this version of
        the rendering which should be a good dead faster.
        res is the width of grid cells

        x - x positions of grid
        y - y positions of grid
        r - particle positions
        rho - particle densities
        bounds - box boundaries xmin,xmax,ymin,ymax
        cutoff - smoothing cutoff

    """

    type = 'lucy'
    nx = x.size
    ny = y.size

    if cutoff > nx: print "cutoff is too large this won't work"

    res = x[1]-x[0] #assume regular
    Z = numpy.zeros((nx,ny))
    xmin,xmax,ymin,ymax = bounds
    #loop over particle positions
    for i in range(m.size):
        # each particle contributes to pixels
        xpixmin = int(floor( (r[i,0] - cutoff - xmin)/res ) )
        ypixmin = int(floor( (r[i,1] - cutoff - ymin)/res ) )
        xpixmax = int(floor( (r[i,0] + cutoff - xmin)/res )   )
        ypixmax = int(floor( (r[i,1] + cutoff - ymin)/res )   )

        # print xpixmin,xpixmax,ypixmin,ypixmax
        # print nx,ny

        for ypx in range(ypixmin,ypixmax):
            for xpx in range(xpixmin,xpixmax):
                
                # where the pixel bounds cross the box bounds, need to 
                # wrap pixel coordinates
                if ypx >= ny: ypx -= ny
                if ypx < 0: ypx += ny
                if xpx >= nx: xpx -= nx
                if xpx < 0: xpx += nx
                
                dr = numpy.array( (r[i] - (x[xpx],y[ypx])))
                # need to apply a minimum image convention
                neighbours.minimum_image(dr,xmax,ymax)
                rsq =  numpy.dot(dr,dr)
                if rsq < (cutoff*cutoff):
                    s = numpy.sqrt(rsq)
                    Z[xpx,ypx] += a[i]*m[i]/rho[i]\
                                * spkernel.lucy_w(s,h[i])
    return Z


def map_to_grid(x,y,r,m,rho,a,h,bounds,cutoff): 
    """ Takes a distribution of sph particles and maps their
        properties onto a regular grid.

        positions (two dimensional array, [particle,spatial])
        property  ([particle])
        two vectors with x and y coordinates
        Outputs
        two dimensional regular array with the smoothed data
        bounds is a tuple (xmin,xmax,ymin,ymax)
    """

    nx = x.size
    ny = y.size
#    X,Y =  pylab.meshgrid(x,y)
    Z = numpy.zeros((nx,ny))
    nlists = []
    ngridpts = x.size*y.size
    xmin,xmax,ymin,ymax = bounds

# instead of building a neighbour list for each pixel, why not build
# a neighbour list as
# (x,y) index
# where (x,y) is the pixel coordinates and index is the neighbour.
# will this be faster?
    
    for i in range(nx):
        for j in range(ny):
            nlist = neighbours.Nlist(ngridpts,cutoff)
            point = numpy.array([x[i],y[j]])
            nlist.build(point,r,h,periodic=True,xmax=xmax,ymax=ymax)
            nlist.assign_properties(m,rho,a,h)
            # compute the sph interpolant
            v = numpy.zeros(len(nlist.ds))
            v = nlist.m/nlist.rho
            # r is the list of distances
            # x is the list of properties
            # vol is the list of volumes
            Z[i,j] = smooth_point(nlist.ds,nlist.dr,nlist.a,v,nlist.h,"lucy")

    return Z 
