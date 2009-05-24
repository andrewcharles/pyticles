"""
    Particles always live in boxes.
    Copyright Andrew Charles 2008
    All rights reserved.
    This code is licensed under the terms of the new BSD license.

"""

XMAX = 64
YMAX = 48

class Box:
    
    def __init__(self,p):
        self.p = p
        self.xmax = XMAX
        self.ymax = YMAX

    def min_image(r):
        print "does not exist yet"

    def apply_periodic_bounds(self,p):
         """ applies periodic boundaries """
         for i in range(p.n):
            if p.r[i,0] > XMAX:
                p.r[i,0] = 0
            if p.r[i,0] < 0:
                p.r[i,0] = XMAX
            if p.r[i,1] > YMAX:
                p.r[i,1] = 0
            if p.r[i,1] < 0:
                p.r[i,1] = YMAX

    def apply_mirror_bounds(self,p):
         """ applies periodic boundaries """
         for i in range(p.n):
            if p.r[i,0] > XMAX:
                p.r[i,0] = XMAX
                p.v[i,0] = -p.v[i,0]
            if p.r[i,0] < 0:
                p.r[i,0] = 0
                p.v[i,0] = -p.v[i,0]
            if p.r[i,1] > YMAX:
                p.r[i,1] = YMAX 
                p.v[i,1] = -p.v[i,1]
            if p.r[i,1] < 0:
                p.r[i,1] = 0
                p.v[i,1] = -p.v[i,1]
