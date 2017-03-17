# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 09:25:03 2017

@author: jarbona
"""
import numpy as np
with open("logp.txt", "w"):
    pass
from collections import namedtuple
SpaceTime = namedtuple("SpaceTime", ["pos", "t"])


class Diffusing:
    def __init__(self, tag):
        self.tag = tag
        self.V = []
        self.bound = []
        self.replicating = []
        self.free = []

    def in_volume_of_interaction(self, ori_list, time):
        for ori in ori_list:
            self.V.append(SpaceTime(ori, time))

    def start_replication(self, ori, time):
        self.replicating.append([SpaceTime(ori, time)])

    def end_replication(self, time):
        self.replicating[-1].append(time)

    def start_bound(self, ori, time):
        self.bound.append([SpaceTime(ori, time)])

    def end_bound(self, time):
        self.bound[-1].append(time)

    def build_time_line(self,maxt=None):
        #Get max_time
        provided = False
        if maxt is None:
            maxt = 0
            for sp in self.V:
                maxt = max(maxt, sp.t)
            for event in self.replicating + self.bound:
                if len(event) == 1:
                    maxt = max(maxt, event[0].t)
                else:
                    maxt = max(maxt, event[1])
        else:
            provided = True

        #print(maxt)
        time_line = np.zeros(int(maxt) + 1)
        for inte in self.V:
            time_line[int(inte.t)] = 1

        for event in self.replicating:
            if len(event) == 1:
                if maxt != event[0].t and not provided:
                    print("Unfinished business")
                    raise
                time_line[int(event[0].t):] = 2
            else:
                time_line[int(event[0].t):int(event[1])] = 2

        for event in self.bound:
            if len(event) == 1:
                if maxt != event[0][1]:
                    print("Unfinished business")
                    raise
                time_line[int(event[0].t):] = 3
            else:
                time_line[int(event[0].t):int(event[1])] = 3
        return time_line




class Origin:

    def __init__(self, tag):
        self.tag = tag
        self.position = self.tag + 0.5
        self.move = False
        self.origin = True
        self.passivated = False
        self.activated = False
        self.path = [[self.tag, 0]]
        self.type = "origin"
        self.activation_time = None

    def state(self):
        add = "O"
        if self.passivated:
            add = "P"
        if self.activated:
            add = "A"
        return add, self.tag, self.path

    def __repr__(self):
        add = ""
        if self.passivated:
            add = "Passi"
        if self.activated:
            add = "Activ"
        return "Origin %i %s" % (self.path[-1][0], add)


class Fork:

    def __init__(self, tag, position, bond_tag, t, diff_diff_tag):
        self.tag = tag  # particle tagg
        self.position = position + 0.5  # Position on the chromosome
        self.bond_tag = bond_tag  # Bond between diffusive elements and monomec
        self.move = True
        self.update_bond = False
        self.origin = False
        self.path = [SpaceTime(int(self.position), t)]
        self.t = t
        self.type = "fork"
        self.diff_diff_tag = diff_diff_tag

    def update_position(self, dt, fork_speed, maxi):

        oldp = int(self.position)

        self.position += self.d * fork_speed * dt
        self.t += dt

        if oldp != int(self.position):

            self.update_bond = True
            # In case of high fork_speed
            while oldp != int(self.position) and oldp >= 0 and oldp <= maxi:
                oldp += self.d
                self.path.append(SpaceTime(oldp, self.t))

        else:
            self.update_bond = False

    def state(self):
        return "F", self.tag, self.path

    def __repr__(self):
        add = "L"
        if self.d == 1:
            add = "R"
        return "%sFork %i" % (add, self.path[-1][0])


class RFork(Fork):

    def __init__(self, tag, position, bond_tag, t, diff_diff_tag):
        Fork.__init__(self, tag, position, bond_tag, t, diff_diff_tag)
        self.d = 1


class LFork(Fork):

    def __init__(self, tag, position, bond_tag, t, diff_diff_tag):
        Fork.__init__(self, tag, position, bond_tag, t, diff_diff_tag)
        self.d = -1


class Polymer():

    def __init__(self, number, start, end, origins):
        self.number = number
        self.start = start
        self.end = end  # End is included
        # print(start, end)
        origins.sort()
        self.origins = origins
        if self.origins != []:
            # print(number,start,np.min(self.origins))
            assert(np.min(self.origins) >= self.start)
            # print(self.end,np.max(self.origins))
            assert(np.max(self.origins) <= self.end)
        self.modules = [Origin(tag) for tag in origins]
        # to keep track of the diff attached in case we attach them one by one
        self.bound_to_origin = {tag: [] for tag in origins}
        # self.replication_state = [0 for i in range(start,end+1)]
        self.t = 0
        self.ended = []
        self.replicated = {r: False for r in range(start, end + 1)}
        self.events = {}

    def has_origin(self, ptag):
        if ptag in self.origins:
            for m in self.modules:
                if m.type == "origin":
                    if m.tag == ptag:
                        return True

            # with open("logp.txt","a") as f:
            #    f.writelines("%i Warning, origin already used %i\n"%(self.number,self.t))
            print("Warning, origin already used 1")
            return False

    def state(self):
        return [self.number, self.start, self.end,
                [m.state for m in self.modules],
                [m.state for m in self.ended]]

    def add_event(self, ptag, event):
        if ptag in self.events:
            self.events[ptag].append([self.t, event])
        else:
            self.events[ptag] = [[self.t, event]]

    def dettach_one_diff(self, ptag, otag):
        if len(self.bound_to_origin[otag]) == 2:
            print("Dettaching one but should have started")
            print(otag, self.bound_to_origin[otag])
            raise
        if self.bound_to_origin[otag][0][0] == ptag:
            self.bound_to_origin[otag] = []
            self.add_event(ptag, "D")

    def attach_one_diff(self, ptag, otag, new_btag):
        if otag not in self.bound_to_origin:
            print("One free origin was not available")
            raise
        else:
            self.bound_to_origin[otag].append([ptag, new_btag])
            self.add_event(ptag, "A")

        if len(self.bound_to_origin[otag]) == 1:
            return False  # Do not start
        else:
            return True

    def get_free_origins(self):
        return self.bound_to_origin.keys()

    def get_n_activ_fork(self):
        return len([m for m in self.modules if not m.origin])

    def get_diff_at_origin(self, otag):
        # list of particles tag , particle bond
        return self.bound_to_origin[otag]

    def add_fork(self, ptags, otag, new_btags, diff_diff_tag):
        found = False
        for i, mod in enumerate(self.modules):
            if mod.origin:
                if mod.tag == otag:
                    found = True
                    break
        if not found:
            print("Warning origin not found")
            raise

        # Not necessary but just for checking elsewhere:
        self.bound_to_origin.pop(otag)
        # with open("logp.txt","a") as f:
        #     f.writelines("%i Warning, origin already used %i\n"%(self.number,self.t))
        self.modules.insert(
            i + 1, RFork(ptags[1], otag, new_btags[1], self.t, diff_diff_tag))
        self.add_event(ptags[1], "R")
        self.modules.insert(
            i, LFork(ptags[0], otag, new_btags[0], self.t, diff_diff_tag))
        self.add_event(ptags[0], "R")

        if self.modules[i + 1].passivated or self.modules[i + 1].activated:
            print("Warning origin already used")
            # with open("logp.txt","a") as f:
            #    f.writelines("%i Warning, origin already used %i\n"%(self.number,self.t))
        self.modules[i + 1].passivated = False
        self.modules[i + 1].activated = True
        self.modules[i + 1].activation_time = self.t
        self.ended.append(self.modules.pop(i + 1))

    def get_replication_profile(self, t=None):
        # self.position_index = range(self.start, self.end + 1)
        self.replication_state = [None for i in range(self.start, self.end + 1)]
        for m in self.modules + self.ended:
            if not m.move:
                continue

            for pos, time in m.path:
                # i = self.position_index.index(pos)
                # print(pos, time)
                self.replication_state[pos - self.start] = time
                # assert(i == pos - self.start)
        return self.replication_state

    def get_interacting_particles(self):
        P = []
        for m in self.ended + self.modules:
            if m.move:
                P.append(m.tag)
        return P

    def get_fork_density(self, fork_speed,cut=0, normed=True,):
        max_t = int(self.t) + int(1 / fork_speed) + 2

        fork_number = np.zeros(max_t)
        rep_p = np.array(self.get_replication_profile())
        for m in self.modules + self.ended:
            if not m.origin:
                # print(m.path)
                if m.path != []:
                    start = m.path[0][1]
                    end = m.path[-1][1]
                    fork_number[int(start) - self.start: int(end) - self.start] += 1

        if normed:
            for t in np.arange(max_t):
                Un_replicated = np.sum(rep_p >= t)
                if Un_replicated == 0:
                    Un_replicated = 1
                fork_number[t] /= Un_replicated
        # print("Done1")
        if cut != 0:
            fork_number[-cut:] = 0
        return fork_number

    def get_norm(self,fork_speed):
        max_t = int(self.t) + int(1 / fork_speed) + 2

        Un_replicated = np.zeros(max_t)

        rep_p = np.array(self.get_replication_profile())
        not_none = ~np.equal(rep_p, None)
        for t in np.arange(max_t):
            # print(not_none)
            # print(rep_p[not_none])
            # print(rep_p)
            Un_replicated[t] = np.sum(rep_p[not_none] >= t) + len(rep_p) - np.sum(not_none)

        return Un_replicated

    def get_DNA_with_time(self, fork_speed):
        # Quantity of replicated dna
        max_t = int(self.t) + int(1 / fork_speed) + 2
        DNA = np.zeros((self.end + 1 - self.start, max_t))
        for m in self.modules + self.ended:
            if not m.move:
                continue
            #if self.number == 0:
        #        print(m.path)
            for pos, time in m.path:
                # i = self.position_index.index(pos)
                # print(pos, time)
                nstep = int(1 / fork_speed) + 2  # To be sure we reach 1 Then fill with 1
                for ctime in range(int(time), int(time) + nstep):
                    #print(ctime)
                    if ctime > int(time):
                        DNA[pos - self.start, min(ctime, max_t - 1)] += min((ctime - int(time)) * fork_speed, 1)
                DNA[pos - self.start, min(ctime, max_t - 1):] = 1

                #print(DNA[pos - self.start])

        #print(DNA[0,::])
        DNA[DNA > 1] = 1
        #print(DNA)
        #if self.number == 0:
    #        print (DNA[:10])
#            print (np.sum(DNA, axis=0)[:10])
        #raise
        return np.sum(DNA, axis=0)

    def get_firing_time_It(self, fork_speed, normed=True, cut=0):

        max_t = int(self.t) + int(1 / fork_speed) + 2

        firing_time = []
        for m in self.modules + self.ended:
            if not m.move:
                if m.activated and m.activation_time is not None:
                    firing_time.append(m.activation_time)
        firing_time.sort()
        It = np.zeros(max_t)
        if cut != 0:
            for el in range(cut):
                firing_time.pop(-1)
        # print(self.t)
        for el in firing_time:
            # print(el,)
            It[int(el)] += 1

        if normed:
            norm = 1.0 * self.get_norm()
            It /= norm
        if cut != 0:
            It[-cut:] = 0
        # print("Done1")
        return firing_time, It

    def get_firing_at_fraction(self, fork_speed, DNA_time, cut=0, bins=100):



        DNA_time /= max(DNA_time)

        firing_time = []
        for m in self.modules + self.ended:
            if not m.move:
                if m.activated and m.activation_time is not None:
                    firing_time.append(DNA_time[int(m.activation_time)])

        firing_time.sort()

        It = np.zeros(bins)

        if cut != 0:
            for el in range(cut):
                firing_time.pop(-1)
        # print(self.t)
        for el in firing_time:
            # print(el,)
            It[min(int(el * bins), bins - 1)] += 1


        return It

    def get_free_origins_time(self, fork_speed, normed=True):
        max_t = int(self.t) + int(1 / fork_speed) + 2

        rep_p = np.array(self.get_replication_profile())
        free = np.zeros(max_t)

        #print("T", int(self.t), self.origins)
        for m in self.modules:
            if not m.move:
                free[:max_t - 1] += 1

        for m in self.ended:
            if not m.move:
                if m.activated:

                    free[:int(m.activation_time)] += 1
                    # print(m.activation_time)
                else:
                    if m.activation_time is None:
                        print("Warning incorrect calculation free origins")
                    else:
                        free[:int(m.activation_time)] += 1
                    # print("Passivated")

                    # print(free)
        if normed:
            for t in np.arange(max_t):
                Un_replicated = np.sum(rep_p >= t)
                if Un_replicated == 0:
                    Un_replicated = 1
                free[t] /= Un_replicated
        # print("Done1")
        return free

        # print(self.modules)

    def increment_time(self, dt, fork_speed, verbose=False):

        self.t += dt
        update_bond = []
        alone = []
        to_release = []
        diff_diff = []
        bind_diff = []
        passivated_origin = []
        if verbose and self.modules != []:
            print(self.start, self.end)
            print(self.modules)
            print(self.ended)
        for m in self.modules:
            if m.move:
                m.update_position(dt=dt, fork_speed=fork_speed, maxi=self.end)
        N_mod = len(self.modules)

        to_remove = []
        #####################################

        im = 0
        while im < N_mod:
            # Take care of passivated Origin Left to right
            m = self.modules[im]
            if m.move:
                if im != N_mod - 1 and m.position >= self.modules[im + 1].position:
                    if self.modules[im + 1].origin:
                        passivated_origin.append(self.modules[im + 1].tag)
                        self.modules[im + 1].passivated = True
                        self.modules[im + 1].activation_time = self.t
                        self.ended.append(self.modules.pop(im + 1))
                        if self.bound_to_origin[self.ended[-1].tag] != []:

                            alone.append(self.bound_to_origin[
                                         self.ended[-1].tag][0][0])
                            if len(self.bound_to_origin[self.ended[-1].tag]) == 2:
                                print("Releasing passivated origin with two diffS")
                        self.bound_to_origin.pop(self.ended[-1].tag)
                        # print("Warning check to release possible attached single diff")
                        # a,b = self.modules[im:im+2]
                        # self.modules[im:im+2] = b,a
                        # im += 1
                        N_mod -= 1

            im += 1

        im = N_mod - 1
        while im > 0:
            # Take care of passivated Origin Right to left
            m = self.modules[im]
            if m.move:
                if im != 0 and m.position <= self.modules[im - 1].position:
                    if self.modules[im - 1].origin:
                        passivated_origin.append(self.modules[im - 1].tag)
                        self.modules[im - 1].passivated = True
                        self.modules[im - 1].activation_time = self.t

                        self.ended.append(self.modules.pop(im - 1))
                        if self.bound_to_origin[self.ended[-1].tag] != []:
                            alone.append(self.bound_to_origin[
                                         self.ended[-1].tag][0][0])
                            # print("Ici")
                            if len(self.bound_to_origin[self.ended[-1].tag]) == 2:
                                print("Releasing passivated origin with two diffS")
                        self.bound_to_origin.pop(self.ended[-1].tag)
                        # print("Warning check to release possible attached single diff")
                        N_mod -= 1
                        # a,b = self.modules[im-1:im+1]
                        # self.modules[im-1:im+1] = b,a
                        # im -= 1

            im -= 1

        for im in to_remove[::-1]:
            m = self.modules.pop(im)
            # print(m.path)

            self.ended.append(m)
            N_mod -= 1
            assert(m.move)

        to_remove = []
        im = 0
        while im < N_mod:
            # Take care of fork collision and bond motion
            m = self.modules[im]
            if m.move:
                if im != N_mod - 1 and m.position >= self.modules[im + 1].position:
                    if self.modules[im + 1].move:
                        # Collision
                        to_release.append(m.bond_tag)
                        to_release.append(self.modules[im + 1].bond_tag)
                        diff_diff.append(m.tag)
                        diff_diff.append(self.modules[im + 1].tag)
                        bind_diff.append([m.tag, self.modules[im + 1].tag])
                        to_remove.append(im)
                        to_remove.append(im + 1)
                        self.add_event(m.tag, "E")
                        self.add_event(self.modules[im + 1].tag, "E")

                        if m.path[-1].pos >= self.modules[im + 1].path[-1].pos:
                            m.path.pop(-1)
                        # if m.path[-1].pos <= self.modules[im + 1].path[-1].pos:
                        #      self.modules[im + 1].path.pop(-1)

                        # Release possible diff_diff bonds
                        if m.diff_diff_tag is not None:
                            to_release.append(m.diff_diff_tag)
                            to_erase = 0 + m.diff_diff_tag
                            # Look for this tag in other fork:
                            for other_fork in self.modules:
                                if other_fork.type == "fork" and other_fork.diff_diff_tag == to_erase:
                                    other_fork.diff_diff_tag = None

                        if self.modules[im + 1].diff_diff_tag is not None:
                            to_release.append(
                                self.modules[im + 1].diff_diff_tag)
                            to_erase = 0 + self.modules[im + 1].diff_diff_tag
                            # Look for this tag in other fork:
                            for other_fork in self.modules:
                                if other_fork.type == "fork" and other_fork.diff_diff_tag == to_erase:
                                    other_fork.diff_diff_tag = None

                        im += 1
                elif m.position < self.start or m.position > self.end:
                    # Take care of fork outside of boundaries
                    alone.append(m.tag)
                    self.add_event(m.tag, "E")
                    # print(alone,m.position)

                    to_release.append(m.bond_tag)

                    # Release possible diff_diff bonds
                    if m.diff_diff_tag is not None:
                        to_release.append(m.diff_diff_tag)
                        to_erase = 0 + m.diff_diff_tag
                        # Look for this tag in other fork:
                        for other_fork in self.modules:
                            if other_fork.type == "fork" and other_fork.diff_diff_tag == to_erase:
                                print("Found")
                                other_fork.diff_diff_tag = None

                    to_remove.append(im)
                    m.update_bond = False  # because we will delete it

                elif m.update_bond:
                    update_bond.append([m.bond_tag, int(m.position)])
                    """    if self.replicated[int(m.position)] and self.number == 0:
                        print("Attaching to already replicated")
                        print(m)
                        print(m.path)
                        print(self.modules)
                        print(self.ended)
                        raise"""
                    self.replicated[int(m.position)] = True
            im += 1

        for im in to_remove[::-1]:
            m = self.modules.pop(im)
            # print(m)
            if m.path != []:
                if m.path[-1][0] < self.start:
                    m.path.pop(-1)
                elif m.path[-1][0] > self.end:
                    m.path.pop(-1)
                # m.path.pop(-1)
                if verbose:
                    print(m.path)
                self.ended.append(m)
            assert(m.move)
        # chek for colisions:
        # for m in
        if verbose and self.modules != []:
            print(self.modules)

        return bind_diff, diff_diff, update_bond, passivated_origin, to_release, alone


if __name__ == "__main__":
    P = Polymer(0, 30, [5, 10, 20])
    P.add_fork([0, 1], 10, ["a", "b"])
    for i in range(11):
        print(P.increment_time(1, verbose=True))
