# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 09:25:03 2017

@author: jarbona
"""

with open("logp.txt","w"):
    pass

class Origin:
    def __init__(self,tag):
        self.tag = tag
        self.position = self.tag
        self.move=False
        self.origin=True
        self.passivated=False
        self.activated=False
        self.path = [[self.tag,0]]
        self.type = "origin"

    def state(self):
        add = "O"
        if self.passivated:
            add="P"
        if self.activated:
            add="A"
        return add,self.tag,self.path
    def __repr__(self):
        add = ""
        if self.passivated:
            add="Passi"
        if self.activated:
            add="Activ"
        return "Origin %i %s"%(self.path[-1][0],add)

class Fork:
    def __init__(self,tag,position,bond_tag,t,diff_diff_tag):
        self.tag = tag  #particle tagg
        self.position = position #Position on the chromosome
        self.bond_tag = bond_tag #Bond between diffusive elements and monomec
        self.move=True
        self.update_bond=False
        self.origin = False
        self.path = [[position,t]]
        self.t = t
        self.type="fork"
        self.diff_diff_tag = diff_diff_tag

    def update_position(self,dt):
        #print(self.position)
        oldp = int(self.position)
        self.position += self.d * dt
        self.t += dt
        if oldp != int(self.position):
            self.update_bond = True
            self.path.append([int(self.position),self.t])
        else:
            self.update_bond =  False

    def state(self):
        return "F",self.tag,self.path

    def __repr__(self):
        add = "L"
        if self.d == 1:
            add="R"
        return "%sFork %i"%(add,self.path[-1][0])

class RFork(Fork):
    def __init__(self,tag,position,bond_tag,t,diff_diff_tag):
        Fork.__init__(self,tag,position,bond_tag,t,diff_diff_tag)
        self.d = 1

class LFork(Fork):
    def __init__(self,tag,position,bond_tag,t,diff_diff_tag):
        Fork.__init__(self,tag,position,bond_tag,t,diff_diff_tag)
        self.d = -1

class Polymer():
    def __init__(self,number,start,end,origins):
        self.number = number
        self.start = start
        self.end = end
        origins.sort()
        self.origins = origins
        self.modules = [Origin(tag) for tag in origins]
        self.bound_to_origin = {tag:[] for tag in origins} # to keep track of the diff attached in case we attach them one by one
        #self.replication_state = [0 for i in range(start,end+1)]
        self.t = 0
        self.ended = []
    def has_origin(self,ptag):
        if ptag in self.origins:
            for m in self.modules:
                if m.type == "origin":
                    if m.tag == ptag:
                        return True


            #with open("logp.txt","a") as f:
            #    f.writelines("%i Warning, origin already used %i\n"%(self.number,self.t))
            print("Warning, origin already used 1")
            return False

    def state(self):
        return [self.number,self.start,self.end,
                [m.state for m in self.modules],
                [m.state for m in self.ended]]

    def attach_one_diff(self,ptag,otag,new_btag):
        if not otag in self.bound_to_origin:
            print("One free origin was not available")
            raise
        else:
            self.bound_to_origin[otag].append([ptag,new_btag])

        if len(self.bound_to_origin[otag]) == 1:
            return False #Do not start
        else:
            return True
    def get_diff_at_origin(self,otag):
        return self.bound_to_origin[otag] # list of particles tag , particle bond

    def add_fork(self,ptags,otag,new_btags,diff_diff_tag):
        found = False
        for i,mod in enumerate(self.modules):
            if mod.tag == otag:
                found = True
                break
        if not found:
            print("Warning origin not found")

        #Not necessary but just for checking elsewhere:
        self.bound_to_origin.pop(otag)
           # with open("logp.txt","a") as f:
           #     f.writelines("%i Warning, origin already used %i\n"%(self.number,self.t))
        self.modules.insert(i+1,RFork(ptags[1],otag,new_btags[1],self.t, diff_diff_tag))
        self.modules.insert(i,LFork(ptags[0],otag,new_btags[0],self.t, diff_diff_tag))
        if self.modules[i + 1].passivated or self.modules[i + 1].activated:
            print("Warning origin already used")
            #with open("logp.txt","a") as f:
            #    f.writelines("%i Warning, origin already used %i\n"%(self.number,self.t))
        self.modules[i+1].passivated=True
        self.modules[i+1].activated=True
        self.ended.append(self.modules.pop(i+1))

    def get_replication_profile(self,t=None):
        self.position_index =  range(self.start,self.end+1)
        self.replication_state = [0 for i in range(self.start,self.end+1)]
        for m in self.modules + self.ended:
            if not m.move:
                continue
            for pos,time in m.path:
                i = self.position_index.index(pos)
                self.replication_state[i] = time
        return self.replication_state

        #print(self.modules)

    def increment_time(self,dt,verbose=False):
        self.t += dt
        update_bond = []
        alone = []
        to_release = []
        diff_diff = []
        bind_diff = []
        passivated_origin = []
        if verbose and self.modules != []:
            print(self.start,self.end)
            print(self.modules)
            print(self.ended)
        for m in self.modules:
            if m.move:
                m.update_position(dt)
        N_mod = len(self.modules)
        im = 0
        to_remove = []
        #####################################
        #Take care of fork outside of boundaries
        while im < N_mod:

            m = self.modules[im]
            if m.move:
                if m.position < self.start or m.position > self.end:
                    alone.append(m.tag)
                    to_release.append(m.bond_tag)

                    #Release possible diff_diff bonds
                    if m.diff_diff_tag is not None:
                        print("La")
                        to_release.append(m.diff_diff_tag)
                        to_erase = 0 + m.diff_diff_tag
                        #Look for this tag in other fork:
                        for other_fork in self.modules:
                            if other_fork.type == "fork" and other_fork.diff_diff_tag == to_erase:
                                print("Found")
                                other_fork.diff_diff_tag = None

                    to_remove.append(im)
                    m.update_bond = False # because we will delete it
                    #Remove last path of the fork if up


            im += 1

        for im in to_remove[::-1]:
            m = self.modules.pop(im)
            m.path.pop(-1)
            self.ended.append(m)
            N_mod -= 1
            assert(m.move)

        to_remove = []
        #####################################

        im = 0
        while im < N_mod:
        #Take care of passivated Origin Left to right
            m = self.modules[im]
            if m.move:
                if im != N_mod -1 and m.position > self.modules[im + 1].position:
                    if self.modules[im + 1].origin:
                        passivated_origin.append(self.modules[im + 1].tag)
                        self.modules[im + 1].passivated = True
                        self.ended.append(self.modules.pop(im + 1))
                        #a,b = self.modules[im:im+2]
                        #self.modules[im:im+2] = b,a
                        #im += 1
                        N_mod -= 1

            im += 1

        im = N_mod-1
        while im > 0:
        #Take care of passivated Origin Right to left
            m = self.modules[im]
            if m.move:
                if im != 0 and m.position < self.modules[im - 1].position:
                    if self.modules[im - 1].origin:
                        passivated_origin.append(self.modules[im - 1].tag)
                        self.modules[im - 1].passivated = True
                        self.ended.append(self.modules.pop(im - 1))
                        N_mod -= 1
                        #a,b = self.modules[im-1:im+1]
                        #self.modules[im-1:im+1] = b,a
                        #im -= 1


            im -= 1


        im = 0
        while im < N_mod:
            #Take care of fork collision and bond motion
            m = self.modules[im]
            if m.move:
                if im != N_mod -1 and m.position > self.modules[im + 1].position:
                    if self.modules[im + 1].move:
                        #Collision
                        to_release.append(m.bond_tag)
                        to_release.append(self.modules[im + 1].bond_tag)
                        diff_diff.append(m.tag)
                        diff_diff.append(self.modules[im + 1].tag)
                        bind_diff.append([m.tag,self.modules[im + 1].tag])
                        to_remove.append(im)
                        to_remove.append(im + 1)

                        #Release possible diff_diff bonds
                        if m.diff_diff_tag is not None:
                            to_release.append(m.diff_diff_tag)
                            to_erase = 0 + m.diff_diff_tag
                            #Look for this tag in other fork:
                            for other_fork in self.modules:
                                if other_fork.type == "fork" and other_fork.diff_diff_tag == to_erase:
                                    other_fork.diff_diff_tag = None

                        if self.modules[im + 1].diff_diff_tag is not None:
                            to_release.append(self.modules[im + 1].diff_diff_tag)
                            to_erase = 0 + self.modules[im + 1].diff_diff_tag
                            #Look for this tag in other fork:
                            for other_fork in self.modules:
                                if other_fork.type == "fork" and other_fork.diff_diff_tag == to_erase:
                                    other_fork.diff_diff_tag = None

                        im += 1
                elif m.update_bond:
                    update_bond.append([m.bond_tag,int(m.position)])
            im += 1
        for im in to_remove[::-1]:
            m = self.modules.pop(im)
            m.path.pop(-1)
            print(m.path)
            self.ended.append(m)
            assert(m.move)
        #chek for colisions:
        #for m in
        if verbose and self.modules != []:
            print(self.modules)
        return bind_diff, diff_diff,update_bond,passivated_origin,to_release,alone
if __name__ == "__main__":
    P = Polymer(0,30,[5,10,20])
    P.add_fork([0,1],10,["a","b"])
    for i in range(11):
        print(P.increment_time(1,verbose=True))
