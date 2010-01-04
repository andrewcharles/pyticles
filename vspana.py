#! /usr/local/bin/python

""" Load a SPAM output file and take a look. (3D)
"""

fname = 'output.nc'
import netCDF4 

ncfile = netCDF4.Dataset(fname,'r')
timev = ncfile.variables['timestep']
r = ncfile.variables['position']
v = ncfile.variables['velocity']
u = ncfile.variables['internal_energy']

# Plot 3D positions of final step
from mpl_toolkits.mplot3d import axes3d
from pylab import figure, show, ion
x = r[-1:,:,:]
x = x.squeeze()
ion()
ax = axes3d.Axes3D(figure())
ax.scatter3D(x[:,0],x[:,1],x[:,2])
show()
print r
