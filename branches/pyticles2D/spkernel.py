#!/usr/bin/python
""" 
    This module provides functions
    that calculate several versions of the
    sph smoothing kernel
    There should be two functions for each kernel. One just gets the kernel
    - this is designed for doing interpolation. The other also gets the kernel gradient - this is designed for doing simulations.
    For example: gauss_kernel, and gauss_kernel_ag (kernel and gradient).
    gauss_kernel_ag will call both the kernel and its derivative.
    Worried a bit about the function call overhead when doing this twice.
"""


from math import *
import numpy

def kernel(r,dx,h,type):
    """ Calculates the kernel and kernel gradient for one pair of particles """
    """ 
    r: distance
    x: displacement vector
    h: smoothing length
    type: kernel type string
    """

    if type == 'gaussian':
        w,dwdx = gauss_kernel(r,dx,h)
    elif type == 'lucy':
        w,dwdx = lucy_kernel(r,dx,h)
    elif type == 'debrun':
        w,dwdx = debrun_kernel(r,dx,h)

    return w,dwdx


def gauss_kernel(r,dx,h):
    """ Ye old Gaussian SPH kernel. Returns w, dwdx """
    """ Works for as many dimensions as you care to inhabit
    """
    
    dx = numpy.atleast_1d(dx)
    dim = dx.size
    q = r/h
    factor = 1.e0 / (h**dim * pi**(dim/2.))
    w = factor * exp(-q*q)
    dwdx=[]
    for i in range(0,dim):
        dwdx[i:] = [w * ( (-2.e0) * (-1.e0) * dx[i]/(h*h))]
    return w,dwdx

def lucy_w(r,h):
    """ The Lucy Kernel. Returns W. Only 2d. """
    q = 5. / (pi * h**2)
    if r < 0:
        r = abs(r)
    if( r < h ):
        w = q * (1 + 3.*r/h)*((1.-r/h))**3
    else:
        w = 0
    return w


def lucy_kernel(r,dx,h):
    """ The Lucy Kernel. Returns W, dwdx"""
    """ Works for one and two dimensions
    """

    dx = numpy.atleast_1d(dx)
    dim = dx.size
    dwdx=[]

    if dim == 1:
        q = 5. / (4. * h)
    elif dim == 2:
        q = 5. / (pi * h**2)

    if r < 0:
        r = abs(r)
    if( r < h ):
        w = q * (1 + 3.*r/h)*((1.-r/h))**3
        if(r == 0):
            dwdx=0.0
        else:
            for i in range(0,dim):
                dwdx[i:] =  [ q * ( (-12./(h**4))*(r**3) +\
                             (24./(h**3))*(r**2) - (12.*r/(h**2)) )* dx[i]/r ]
    else:
        w = 0
        for i in range(0,len(dx)):
            dwdx[i:] = [0]
    return w,dwdx

def debrun_kernel(r,dx,h):
    """ The Debrun spiky kernel. Returns a tuple of w, dwdx """
    """ Works for two dimensions only
        Is not normalised correctly for one dimension.
        (in a big way)
    """

    #! gamma_g := 20/(Pi*h^6);
    #! w_spiky_2d := (r,h) -> piecewise((r<(h)), gamma_g * (h-r)^3, 0);
    dx = numpy.atleast_1d(dx)
    dim = dx.size
    dwdx=[] 
    if r < 0:
        r = abs(r)
    if( r < h ):
        if dim == 1:
            q = ( ( (-7.*(h**4)/4) + 2*(h**4) )**(-1) )/ 2  
            #(2.0/(3.0*(pi*(h**2)))) 
        elif dim ==2:
            q = (15.0/(pi*(h**6)))
        w =  q * ((h - r)**3)
        if(r == 0):
            dwdx=0.0
        else:
            for i in range(0,dim):
                dwdx[i:] = [q * (-3.0*(h**2) + 6.0*h*r - 3.0*(r**2)) \
                              * (dx[i]/r) ]
    else:
        for i in range(0,dim):
            dwdx[i:] = [0]
        w = 0
    return w,dwdx

def test():
    """ Tests the module """
    print "Testing spkernel.py"

    print "Testing for one dimension"
    r = 1.0
    dx = 1.0
    h = 2.0
    type = "gaussian"
    
    print "Gaussian"
    print kernel(r,dx,h,"gaussian")
    print gauss_kernel(r,dx,h)
    
    print "Debrun"
    print kernel(r,dx,h,"debrun")
    print debrun_kernel(r,dx,h)
    
    print "Lucy"
    print kernel(r,dx,h,"lucy")
    print lucy_kernel(r,dx,h)

    print "Testing for two dimensions"
    r = 1.4
    dx = 1.0,1.0
    h = 2.0
    
    print "Gaussian"
    print kernel(r,dx,h,"gaussian")
    print gauss_kernel(r,dx,h)

    print "Debrun"
    print kernel(r,dx,h,"debrun")
    print debrun_kernel(r,dx,h)
    
    print "Lucy"
    print kernel(r,dx,h,"lucy")
    print lucy_kernel(r,dx,h)

    # Test for three dimensions

if __name__ == "__main__":
    test()