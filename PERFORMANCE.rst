A place for notes on making this version very fast.

forcetime - time for all forces (Called in derivatives)
derivtime - time for the first derivatives call before calling integrate
pairtime - time for pairseps for both neighbour lists (Called in derivatives)
integtime -  time for the integration step. This calls derivatives
spamtime - time for one call of spam properties (computes for two neighbour lists, but not all properties for both.) (Called in derivatives)
updatetime - Full update time (dominated by integtime)

integtime (for RK4) should be 4 * derivtime

derivtime = forcetime + spamtime + pairtime

So it's forcetime, spamtime and pairtime that need to be examined.


20091116
For 27 particles, 359 pairs
forcetime - 0.016
derivtime - 0.028
pairtime - 0.011
integtime - 0.113
spamtime - 0.001

By these figures, forcetime and pairtime are the weak links in the chain.
Fortranifying the neighbour list might help, but this would make the interface
considerably more complex. Algorithmic improvements to the list may be
a better use of resources.

By using c_forces.pyx

forcetime - 0.000

Wow. Now i really need to hit the neighbour list. I've just made a new pairsep
pyx module.

pairtime - 0.004

Not that big an increase - perhaps there are things I'm not getting right
I changed the loop variables to cython unsigned ints. This led to

pairtime - 0.002

about an order of magnitude reduction, but still not spectacular. Maybe it's the object references?

2010-August

Still it's too slow for sytems of thousands of particles.
The most time now is spent on force computation, although pair seperation is
also major, at about half the time of force compute. Weighing up whether to write
a fortran front end (sp3d) to get the optimisation going, or work with python.

One place to look at is the wraptest.py script (consider changing the name to
something more descriptive... )




