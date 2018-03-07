# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 16:39:18 2017

@author: jarbona
"""


import numpy as np
from cpolymer.polymer import Polymer
from cpolymer.lsimu import LSimu
from cpolymer.constrain import Sphere, Point

from cpolymer.halley.constrain import Spherical
from cpolymer.halley.vectors import V

len_chrom = [46, 162, 63, 306, 115, 54, 218, 112, 87, 149, 133, 215, 184, 156, 218, 189]
dist_centro = [30, 47, 22, 89, 30, 29, 99, 21, 71, 87, 88, 30, 53, 125, 65, 111]
p_ribo = [[0, 0] for i in range(16)]
p_ribo[11] = [90, 150]

Radius = 16.6
Mt = 0.4 * 16.6


dnuc = 3
bead_type = 1  # All the beads are going to be the same type
liaison = {"1-1": [1, 1], "1-2": [1, 2], "1-3": [1, 3], "1-4": [(dnuc + 1) / 2., 4], "1-5": [0, 5],
           "2-2": [1, 6], "2-3": [1, 7], "2-4": [(dnuc + 1) / 2., 8], "2-5": [0, 9],
           "3-3": [1, 10], "3-4": [(dnuc + 1) / 2., 11], "3-5": [Mt, 12],
           "4-4": [dnuc, 13], "4-5": [0, 14],
           "5-5": [0, 15]}


def create_init_conf_yeast(len_chrom=len_chrom, p_ribo=p_ribo, dist_centro=dist_centro, Radius=Radius, Mt=Mt, dnuc=dnuc):

    Sim = LSimu()

    Nchromosomes = len(len_chrom)
    nucleus = Sphere(position=[0, 0, 0], radius=Radius)

    for X in range(Nchromosomes):
        # print(X)
        # We need to define geometrical constrain:
        # The centromere must be at a distance mt of the spb positioned at (-Radius,0,0)
        # We use the module halley for that
        S_spb = Spherical(V(-Radius, 0, 0), radius=Mt)  # we define a sphere centered on the spb
        # a sphere centered that intersect the spb sphere
        S_centered = Spherical(V(0, 0, 0), radius=Radius - Mt * 0.85)

        circle = S_spb * S_centered
        centromere = circle.get_random()

        # We must then construct a sphere centered on the centromere with a radius sqrt(Nbead)
        # and look at its intersection with the nucleus
        # a sphere centered that intersect the spb sphere
        Nucleus = Spherical(V(0, 0, 0), radius=Radius * 0.95)
        d1 = dist_centro[X]
        Telo1_possible = Spherical(centromere, radius=np.sqrt(d1)) * Nucleus
        inc = 0
        while Telo1_possible.__class__.__name__ == "Nowhere":
            inc += 1
            Telo1_possible = Spherical(centromere, radius=np.sqrt(d1 + inc)) * Nucleus
        telo1 = Telo1_possible.get_random()

        d2 = len_chrom[X] - dist_centro[X]
        Telo2_possible = Spherical(centromere, radius=np.sqrt(d2)) * Nucleus
        inc = 0
        while Telo2_possible.__class__.__name__ == "Nowhere":
            inc += 1
            Telo2_possible = Spherical(centromere, radius=np.sqrt(d2 + inc)) * Nucleus
        telo2 = Telo2_possible.get_random()

        if X != 11 or (X == 11 and p_ribo[11][1] == 0):
            Sim.add(Polymer(N=len_chrom[X], type_bead=[2] + [1] * (d1 - 2) + [3] + [1] * (d2 - 1) + [2],
                            liaison=liaison,
                            gconstrain=[nucleus],
                            rigid_constrain=False,
                            lconstrain=[Point(index=0, position=telo1._v),
                                        Point(index=d1, position=centromere._v),
                                        Point(index=len_chrom[X] - 1, position=telo2._v)],
                            rc=0.5))
        else:
            # This chromosome is different because it has a nucleole
            p_r, size = p_ribo[11]
            Beads = [2] + [1] * (d1 - 2) + [3] + [1] * (p_r - d2) + [4] * \
                size + [1] * (len_chrom[X] - p_r - 1) + [2]
            # print(Beads)
            print(len(Beads),len_chrom[X] + size,d1,d2)
            Sim.add(Polymer(N=len_chrom[X] + size, type_bead=Beads,
                            liaison=liaison,
                            rigid_constrain=False,
                            gconstrain=[nucleus],
                            lconstrain=[Point(index=0, position=telo1._v),
                                        Point(index=d1, position=centromere._v),
                                        # We add a new constrain: the center
                                        Point(index=int(p_r + size / 2),
                                              position=(0.66 * Radius, 0, 0)),
                                        # of the nucleole must be at 2/3 of
                                        # the radius at the
                                        # opposite of
                                        # the spb:
                                        Point(index=len_chrom[X] - 1, position=telo2._v)],
                            rc=0.5))

    return Sim
