""" Time the nlist creation and pair seperation.
    Create random particle positions for 1000 particles
    plot_nlist_time.py will plot the data generated by
    this script.
"""
import particles
import neighbour_list
from time import time
import numpy as np

for j in range(1):

        fname = 'pairtime' + str(j) + '.dat'
        ofile = open(fname,'w')

        for n in range(3,100):
                p = particles.ParticleSystem(n,d=3,maxn=n)
                for i in range(n):
                        p.r[i,:] = np.random.random(3) * 10
                nl = neighbour_list.VerletList(p,cutoff=10,tolerance=2)
                nl.build()
                nl.compress()
                t = time()
                nl.separations()
                # execution time is in seconds
                ofile.write('%d %d %f\n' %(n,nl.nip,(time() - t)) )

        ofile.close()



