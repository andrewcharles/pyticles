""" Generating configurations of particle positions.
   
    randpt --
    grid --
    grid3d --
    random -- 

"""
import numpy

def randpt(xmin,xmax,ymin,ymax):
    """ returns a numpy array representing a random point """
    x = numpy.random.rand() * (xmax - xmin) + xmin
    y = numpy.random.rand() * (ymax - ymin) + ymin
    return numpy.array((x,y))

def grid(n,xside,yside,origin,spacing=1.0):
    """ Makes a regular spaced grid with the sides
        as specified.
        Origin is a tuple x,y coordinate.
    """
    x = numpy.zeros([n,2])
    j,k = origin
    for i in range(n):
        x[i,0]=j
        x[i,1]=k
        j=j+spacing
        if (i+1)%xside == 0:
            j=origin[0]
            k=k+spacing
    return x

def grid3d(n,side,origin,spacing=1.0):
    """ Makes a regular spaced grid with the sides
        as specified.
        side -- x,y,z tuple: number of points along each side
        origin -- a tuple x,y,z coordinate.
    """
    r = numpy.zeros([n,3])
    x,y,z = origin
    q = 0
    for i in range(side[0]):
        for j in range(side[1]):
            for k in range(side[2]):
                r[q,0] = (x - (side[0]-1)*spacing/2.) + i * spacing
                r[q,1] = (y - (side[1]-1)*spacing/2.) + j * spacing
                r[q,2] = (z - (side[2]-1)*spacing/2.) + k * spacing
                q += 1
                if q >= n:
                    return r
    return r

def hotspotgrid3d(n,side,origin,spacing=1.0,temp=(0.2,1.9)):
    """ Makes a regular spaced grid with the sides
        as specified. The central third of points are hot.
        side -- x,y,z tuple: number of points along each side
        origin -- a tuple x,y,z coordinate.
    """
    r = numpy.zeros([n,3])
    t = numpy.zeros([n])
    x,y,z = origin
    q = 0
    # get the thirds
    xt = side[0]/3,2*side[0]/3
    yt = side[1]/3,2*side[1]/3
    zt = side[2]/3,2*side[2]/3
    if side[2] == 1: zt = 0,0
    print xt,yt,zt
    for i in range(side[0]):
        for j in range(side[1]):
            for k in range(side[2]):
                r[q,0] = (x - (side[0]-1)*spacing/2.) + i * spacing
                r[q,1] = (y - (side[1]-1)*spacing/2.) + j * spacing
                r[q,2] = (z - (side[2]-1)*spacing/2.) + k * spacing
                if i>=xt[0] and i<=xt[1] and j>=yt[0] and j<=yt[1] \
                    and k>=zt[0] and k<=zt[1]:
                    t[q] = temp[1]
                else:
                    t[q] = temp[0]
                q += 1
                if q >= n:
                    return r,t
    return r,t
    
def random(n,xmin,xmax,ymin,ymax):
    """ 
    """
    r = numpy.zeros((n,2))
    width = xmax - xmin
    height = ymax - ymin
    for i in range(n):
       r[i,:] = numpy.random.random()*width + xmin,\
                numpy.random.random()*height + ymin 
    return r


def test():
    from mpl_toolkits.mplot3d import axes3d
    from pylab import figure, show, ion
    r = grid3d(1000,(10,10,10),(0.5,0.5,0.5),spacing=0.1)
    ax = axes3d.Axes3D(figure())
    ax.scatter3D(r[:,0],r[:,1],r[:,2])
    ion()
    show()
    print r


if __name__ == '__main__':
    test()

