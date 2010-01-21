#! /usr/local/bin/python

""" Load a SPAM output file and take a look. (3D)
"""

fname = 'output.nc'
#fname = 'op.nc'
import netCDF4 
import numpy as np
from time import sleep

ncfile = netCDF4.Dataset(fname,'r')
timev = ncfile.variables['timestep']
r = ncfile.variables['position']
v = ncfile.variables['velocity']
u = ncfile.variables['internal_energy']

# Plot 3D positions of final step
from mpl_toolkits.mplot3d import axes3d
from pylab import figure, show, ion, clf, gcf
ion()
ax = axes3d.Axes3D(figure())


def plotstep(step):
    x = r[step,:,:]
    if np.isnan(x).any():
        print 'NaN present'
    else:
        x = x.squeeze()
        ax.scatter3D(x[:,0],x[:,1],x[:,2])
        #show()

def clear():
    global ax
    clf()
    ax = axes3d.Axes3D(gcf())

def plotall():
    for i in range(10):
        x = r[i,:,:]
        if np.isnan(x).any():
            print 'NaN present'
            break
        else:
            x = x.squeeze()
            ax.scatter3D(x[:,0],x[:,1],x[:,2])
            sleep(2.0)
