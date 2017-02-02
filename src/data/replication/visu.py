# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 12:42:19 2017

@author: jarbona
"""



from mayavi import mlab
from tvtk.tools import visual
# from colour import Color
import mdtraj as md
import numpy as np
try:
    import _pickle as cPickle
except:
    import cPickle
import copy
import sys
import os
# Create a figure
f = mlab.figure(size=(500,500))
# Tell visual to use this as the viewer.
visual.set_viewer(f)

# A silly visualization.
color_d ={}
color_d["green"] = (0,128,0)
# Even sillier animation.

data_folder = sys.argv[1]

t = md.load(os.path.join(data_folder, 'poly.dcd'),
            top=os.path.join(data_folder, "atoms.hoomdxml"))

class Representation:
    def __init__(self,traj,what,by_resid=True,tube=True,t=0,update=False,time_color=None,color="red"):
        self.first = True
        self.traj = traj
        self.selections = []
        self.time_color = time_color
        self.tube = tube
        self.color=color
        self.named_selections = []

        n_chains = len(list(self.traj.topology.chains))


        if by_resid:
            for i in range(n_chains):
                add = ""
                if what != "":
                    add = " and %s"%what
                sele = "resid %i"%i+add

                self.named_selections.append(sele)
                sel = self.traj.topology.select(sele)

                #print(sele,len(sel))
                if len(sel) != 0:
                    self.selections.append(sel)
        else:
            self.named_selections.append(what)
            sel = self.traj.topology.select(what)
            if len(sel) != 0:
                self.selections.append(sel)
            else:
                print("Empty selection",what)

    def draw(self,time):
        #print(self.selections)
        if self.first:

            self.rep = []

        for isel,sel in enumerate(self.selections):
            x,y,z = self.traj.xyz[time, sel].T
            if self.time_color != None:


                colors = self.time_color(time)
                if len(colors) < isel + 1:
                    print("The colouring function did not cover all the selection")
                    print(isel,self.named_selections[isel])
                    print(len(x),isel)
                    raise

                colors = colors[isel]

                if len(x) != len(colors):
                    print(len(x),len(colors),isel)
                    print(self.named_selections[isel])
                    raise
            else:
                colors = [1 for i in range(len(x))]

            #print(len(x),len(colors))

            if self.first:
                if self.tube:
                    self.rep.append(mlab.plot3d(x, y, z,colors ,vmin=1,vmax=2.1,tube_radius=0.03,opacity=1))
                else:
                    self.rep.append(mlab.points3d(x, y, z,colors ,scale_factor=0.1,color=color_d.get(self.color, "green")))
            else:

                self.rep[isel].mlab_source.set(x=x,y=y,z=z,scalars=colors)#,scalars=replic)
        self.first = False
#print(x.shape,y.shape)

with open(os.path.join(data_folder,"polymer_timing.dat"),"rb") as f:
        lPolymers = cPickle.load(f)


P = [np.array(lPolymers[p].get_replication_profile()) for p in range(len(lPolymers))]

def time_color(time):
    replic = copy.deepcopy(P)
    for l in range(len(P)):
        replic[l][replic[l] > time / 10. ] = 0
        replic[l][replic[l] != 0 ] = 2
        replic[l][replic[l] == 0 ] = 1
    #print(len(replic))
    return replic

Reprs = [Representation(traj=t,what="name Diff",tube=False,by_resid=False),
         Representation(traj=t,what="name Telo",tube=False,by_resid=False,color="green"),
         Representation(traj=t,what=" (not name Diff) and (not name Spb)",time_color=time_color)]
for r in Reprs:
    r.draw(0)

#Load the replication information:


i = 0
@mlab.show
@mlab.animate(delay=100)
def anim():
    """Animate the b1 box."""
    i = 0
    while 1:

        for r in Reprs:
            r.draw(i)
        i = i+ 1
        i = i % len(t)
        yield

# Run the animation.
anim()
