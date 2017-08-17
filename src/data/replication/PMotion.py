# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 09:25:03 2017

@author: jarbona
"""
import numpy as np
with open("logp.txt", "w"):
    pass
from collections import namedtuple
SpaceTime = namedtuple("SpaceTime", ["ori", "t", "pos"])


class SpaceTimeB:

    def __init__(self, pos, time, outpos, bead):
        self.pos = pos
        self.time = time
        self.outpos = outpos
        self.bead = bead

    def get_values(self):
        return self.pos, self.time, self.outpos, self.bead

    def __repr__(self):
        if self.outpos is not None:
            return "(t=%.2f, p=%.1f, otp=%.1f, bead= %.1f)" % (self.time, self.pos, self.outpos, self.bead)
        else:
            return "(t=%.2f, p=%.1f)" % (self.time, self.pos)


class Diffusing:

    def __init__(self, tag):
        self.tag = tag
        self.V = []
        self.bound = []
        self.replicating = []
        self.free = []
        self.isfree = True
        self.isreplicating = False
        self.isbound = False

    def change_state(self, name):
        names = ["isfree", "isreplicating", "isbound"]
        for inames in names:
            setattr(self, inames, False)
        setattr(self, "is" + name, True)

    def in_volume_of_interaction(self, ori_list, time):
        for ori in ori_list:
            self.V.append(SpaceTime(ori, time))

    def start_replication(self, ori, time, pos=None):
        if hasattr(self, "change_state"):
            self.change_state("replicating")
        self.replicating.append([SpaceTime(ori, time, pos=pos)])

    def end_replication(self, time, pos=None):
        if hasattr(self, "change_state"):
            self.change_state("free")
        self.replicating[-1].append(SpaceTime(self.replicating[-1][0].ori, time, pos=pos))

    def start_bound(self, ori, time, pos=None):
        if hasattr(self, "change_state"):
            self.change_state("bound")
        self.bound.append([SpaceTime(ori, time, pos=pos)])

    def end_bound(self, time, pos=None):
        if hasattr(self, "change_state"):
            self.change_state("free")
        self.bound[-1].append(SpaceTime(self.bound[-1][0].ori, time, pos=pos))

    def distance_between_boundings(self, time=False, only_replicating=False, pos=None, real_d=True):
        # sort events then get distances
        if only_replicating:
            merged = sorted(self.replicating, key=lambda x: x[0].t)
        else:
            merged = sorted(self.bound + self.replicating, key=lambda x: x[0].t)
        dists = []
        if len(merged) >= 2:
            for event2, event1 in zip(merged[1:], merged[:-1]):
                if time:
                    dists.append(event2[0].t - event1[1].t)
                else:
                    if pos is None:
                        if real_d:
                            dists.append(np.linalg.norm(
                                np.array(event2[0].pos) - np.array(event1[1].pos)))
                        else:
                            dists.append(np.linalg.norm(
                                np.array(event2[0].pos) - np.array(event1[0].pos)))
                    else:
                        dists.append(np.linalg.norm(
                            np.array(pos[event2[0].ori]) - np.array(pos[event1[0].ori])))

        return np.array(dists)

    def pos_boundings(self, time=False, only_replicating=False):
        # sort events then get distances
        if only_replicating:
            merged = sorted(self.replicating, key=lambda x: x[0].t)
        else:
            merged = sorted(self.bound + self.replicating, key=lambda x: x[0].t)

        pos = []
        if len(merged) >= 1:
            for event in merged:
                if time:
                    pos.append(event[0].t)
                    pos.append(event[1].t)
                else:
                    # print(event)
                    pos.append(event[0].pos)
                    pos.append(event[1].pos)
        return np.array(pos)

    def build_time_line(self, maxt=None):
        # Get max_time
        provided = False
        if maxt is None:
            maxt = 0
            for sp in self.V:
                maxt = max(maxt, sp.t)
            for event in self.replicating + self.bound:
                if len(event) == 1:
                    maxt = max(maxt, event[0].t)
                else:
                    maxt = max(maxt, event[1].t)
        else:
            provided = True

        # print(maxt)
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
                time_line[int(event[0].t):int(event[1].t)] = 2

        for event in self.bound:
            if len(event) == 1:
                if maxt != event[0].t:
                    print("Unfinished business")
                    raise
                time_line[int(event[0].t):] = 3
            else:
                time_line[int(event[0].t):int(event[1].t)] = 3
        return time_line


class Origin:

    def __init__(self, tag, random=False, position=None, strength=1):
        self.tag = tag
        if position is not None:
            self.position = position
        else:
            if random:
                self.position = self.tag + 0.0001 + 0.9997 * np.random.rand()
            else:
                self.position = self.tag + 0.5
        # print(tag,self.position,position)
        self.move = False
        self.origin = True
        self.passivated = False
        self.activated = False
        self.path = [[self.position, 0]]
        self.type = "origin"
        self.activation_time = None
        self.strength = strength

    def state(self):
        add = "O"
        if self.passivated:
            add = "P"
        if self.activated:
            add = "A"
        return add, self.tag, self.path

    def get_strength(self):
        return self.strength

    def __repr__(self):
        add = ""
        if self.passivated:
            add = "Passi"
        if self.activated:
            add = "Activ"
        return "Origin %.1f %s" % (self.position, add)


def cut_path(start, end, direction):
    initd = 0 + direction
    if direction < 0:
        start, end = end, start
        direction = -1 * direction

    initpos = 0 + start
    delta = end - start
    path = [0 + initpos]
    cond = lambda x: x <= end

    while (initpos + delta) != int(initpos) and cond(initpos):
        ddelta = int(initpos) + direction - initpos
        initpos += ddelta
        ddelta -= initpos
        path.append(initpos)
    path[-1] = end
    if len(path) >= 2 and path[-1] == path[-2]:
        path.pop(-1)

    if initd < 0:
        path = path[::-1]
    return path


class Fork:

    def __init__(self, tag, position, bond_tag, t, diff_diff_tag):
        self.tag = tag  # particle tagg
        self.position = position  # Position on the chromosome
        self.bond_tag = bond_tag  # Bond between diffusive elements and monomec
        self.move = True
        self.update_bond = False
        self.origin = False
        self.path = [SpaceTimeB(self.position, t, None, int(self.position))]
        self.t = t
        self.type = "fork"
        self.diff_diff_tag = diff_diff_tag

    def update_position(self, dt, fork_speed, mini, maxi):

        # print("Pos",self.position,self.__repr__())
        old_realp = 0 + self.position
        oldp = int(self.position)
        self.position += self.d * fork_speed * dt
        self.t += dt
        if (self.d > 0 and oldp != int(self.position)) or (self.d < 0 and (oldp != int(self.position) or self.position == self.path[-1].bead)):
                    # if oldp != int(self.position) or
            self.update_bond = True
            path = cut_path(old_realp, self.position, self.d)
            init_t = 0 + self.t - dt

            if self.d > 0:
                if int(path[-1]) != path[-1] > 0:
                    path.pop(-1)
                for ip, p in enumerate(path[1:]):
                    self.path[-1].outpos = p
                    init_t += abs(p - path[ip]) / fork_speed
                    if int(p + 0.1 * self.d) <= maxi + 1:
                        self.path.append(SpaceTimeB(p, init_t, None, int(p + 0.1 * self.d)))
                    else:
                        break

            else:
                if int(path[-1]) != path[-1] > 0:
                    path.pop(-1)
                for ip, p in enumerate(path[1:]):
                    self.path[-1].outpos = p
                    init_t += abs(p - path[ip]) / fork_speed
                    if int(p + 0.1 * self.d) >= mini and int(p + 0.1 * self.d) != self.path[-1].bead:
                        self.path.append(SpaceTimeB(p, init_t, None, int(p + 0.1 * self.d)))
                    else:
                        break

            # Find time of entrance:
            """
            if self.d == 1:


                entrance = int(self.position)
                time = self.t - (self.position - entrance) / fork_speed
                bead = entrance

            else:
                entrance = int(self.position) + 1
                time = self.t + (self.position - entrance) / fork_speed
                bead = entrance - 1
            # In case of high fork_speed
            # while oldp != int(self.position) and oldp >= 0 and oldp <= maxi:
            #    oldp += self.d
            self.path[-1].outpos = entrance
            self.path.append(SpaceTimeB(entrance,time,None,bead))"""

        else:
            self.update_bond = False

    def state(self):
        return "F", self.tag, self.path

    def __repr__(self):
        add = "L"
        if self.d == 1:
            add = "R"
        return "%sFork %.1f" % (add, self.position)


class RFork(Fork):

    def __init__(self, tag, position, bond_tag, t, diff_diff_tag):
        Fork.__init__(self, tag, position, bond_tag, t, diff_diff_tag)
        self.d = 1


class LFork(Fork):

    def __init__(self, tag, position, bond_tag, t, diff_diff_tag):
        Fork.__init__(self, tag, position, bond_tag, t, diff_diff_tag)
        self.d = -1


class Polymer():

    def __init__(self, number, start, end, origins, random=False, positions=None, strengths=None):
        self.number = number
        self.start = start
        self.end = end  # End is included
        # print(start, end)
        origins.sort()
        self.origins = origins
        if self.origins != []:
            # print(number,start,np.min(self.origins))
            # assert(np.min(self.origins) >= self.start)
            # print(self.end,np.max(self.origins))
            # assert(np.max(self.origins) <= self.end)
            pass

        if strengths is None:
            strengths = np.ones_like(self.origins)

        if positions is not None:
            if len(positions) != len(origins):

                print("Error on position,PMotion.py", len(positions), len(origins))
                raise
            self.modules = [Origin(tag, random=random, position=position, strength=strength)
                            for tag, position, strength in zip(origins, positions, strengths)]
        else:
            self.modules = [Origin(tag, random=random, strength=strength)
                            for tag, strength in zip(origins, strengths)]
        #print(origins, strengths)
        self.o_strength = {tag: strength for tag, strength in zip(origins, strengths)}
        # to keep track of the diff attached in case we attach them one by one
        self.bound_to_origin = {tag: [] for tag in origins}
        self.DNA_state = np.arange(self.start, self.end + 1)
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
            print("Warning origin not found", otag)
            raise

        # Not necessary but just for checking elsewhere:
        self.bound_to_origin.pop(otag)
        # with open("logp.txt","a") as f:
        #     f.writelines("%i Warning, origin already used %i\n"%(self.number,self.t))
        self.modules.insert(
            i + 1, RFork(ptags[1], mod.position, new_btags[1], self.t, diff_diff_tag))
        self.add_event(ptags[1], "R")
        self.modules.insert(
            i, LFork(ptags[0], mod.position, new_btags[0], self.t, diff_diff_tag))
        self.add_event(ptags[0], "R")

        if self.modules[i + 1].passivated or self.modules[i + 1].activated:
            print("Warning origin already used")
            # with open("logp.txt","a") as f:
            #    f.writelines("%i Warning, origin already used %i\n"%(self.number,self.t))
        self.modules[i + 1].passivated = False
        self.modules[i + 1].activated = True
        self.modules[i + 1].activation_time = self.t
        self.ended.append(self.modules.pop(i + 1))

    def get_replication_profile(self, fork_speed, t=None, dt=1):
        prof = self.get_DNA_with_time(fork_speed=fork_speed, dt=dt)[1]
        # print("la")
        return np.argmax(prof, axis=1)
        # self.position_index = range(self.start, self.end + 1)

    def get_interacting_particles(self):
        P = []
        for m in self.ended + self.modules:
            if m.move:
                P.append(m.tag)
        return P

    def get_fork_density(self, fork_speed, cut=0, normed=False, dt=1):
        max_t = int(self.t / dt) + int(1 / fork_speed / dt) + 2

        fork_number = np.zeros(max_t)
        # rep_p = np.array(self.get_replication_profile(fork_speed=fork_speed))
        for m in self.modules + self.ended:
            if not m.origin:
                # print(m.path)
                if m.path != []:
                    start = int(round(m.path[0].time / dt, 0))
                    end = int(round(m.path[-1].time / dt, 0))
                    fork_number[int(start): int(end)] += 1

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

    def get_dist_between_activated_origins(self, fork_speed, time=None, dt=1):

        max_t = int(self.t / dt) + int(1 / fork_speed / dt) + 2

        firing_time = []
        for m in self.modules + self.ended:
            if not m.move:
                if m.activated and m.activation_time is not None:
                    if time is not None and m.activation_time / 1.0 / dt > time:
                        continue
                    firing_time.append([int(round(m.activation_time / 1.0 / dt, 0)), m.position])
        firing_time.sort(key=lambda x: x[1])

        firing_position = np.array(firing_time)

        Dist = []
        for time in range(max_t):
            fired = firing_position[::, 0] <= time
            dist = firing_position[fired][::, 1]
            dist = dist[1:] - dist[:-1]
            Dist.append(dist)
        return Dist, firing_position

    def get_norm(self, fork_speed):
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

    def get_correlations(self, fork_speed, dt=1, thresh=1, merge=0):

        DNA_time = self.get_DNA_with_time(fork_speed=fork_speed, dt=dt)[1]
        # print(DNA_time.shape)
        # axis0 = length axis1 = time

        def return_connected_slow(a):
            # print(a)
            l = [[]]
            for i, el in enumerate(a):
                if (el is False or el == 0) and l[-1] != []:
                    l.append([])
                if (el is True or el == 1):
                    l[-1].append(i)
            if l[-1] == []:
                l.pop(-1)
            return l

        def return_connected(a):
            boo = np.array(a, dtype=np.bool)
            indices = np.nonzero(boo[1:] != boo[:-1])[0] + 1
            # print(indices)
            b = np.split(np.arange(len(a)), indices)
            b = b[0::2] if boo[0] else b[1::2]
            return b

        def IOD_IRTD_TL(a, thresh=1):
            # IOD = inter origin distance
            # IRTD = inter replication track distance
            # TL = track length
            toa = a >= thresh

            connected = return_connected(toa)
            # print(connected)

            TL = [len(component) for component in connected]
            IRTD = []
            IOD = []
            if len(connected) > 2:
                for c1, c2 in zip(connected[:-1], connected[1:]):
                    IRTD.append(c2[0] - c1[-1])

                CM = np.array([np.mean(component) for component in connected])
                # for c1, c2 in zip(CM[:-1], CM[1:]):
                #    IOD.append(c2 - c1)
                IOD = CM[1:] - CM[:-1]
            return IOD, IRTD, TL

        IODs = []
        IRTDs = []
        TLs = []

        for struct in DNA_time.T:
            # print(struct.shape)
            # print(struct)
            iod, irtd, tl = IOD_IRTD_TL(struct, thresh=thresh)
            IODs.append(iod)
            IRTDs.append(irtd)
            TLs.append(tl)
        return IODs, IRTDs, TLs

    def get_DNA_with_time(self, fork_speed, dt=1):
        # Quantity of replicated dna
        if hasattr(self, "cache"):
            if int(self.t / dt) + int(1 / fork_speed / dt) + 2 == self.max_t:
                return self.sDNA, self.DNA
        max_t = int(self.t / dt) + int(1 / fork_speed / dt) + 2
        DNA = np.zeros((self.end + 1 - self.start, max_t))
        for m in self.modules + self.ended:
            if not m.move:
                continue
            # if self.number == 0:
        #        print(m.path)

            # pos = m.path[0].pos - m.d
        #    print(m.path)
            # print(m.path)

            for ip, spaceTimeB in enumerate(m.path):
                pos, time, outpos, bead = spaceTimeB.get_values()

                if ip > 1 and pos != m.path[ip - 1].outpos:
                    print(m.path)

                # print(pos,outpos,time)

                outtime = (time + abs(outpos - pos) / fork_speed)

                time /= dt
                outtime /= dt
                # print(outtime)
                times = [time, outtime]
                t = int(round(time, 0)) + 1
                # print
                while t < times[-1]:
                    times.insert(-1, t)
                    t += 1
                # print("t",times)

                for t, tp1 in zip(times[:-1], times[1:]):
                    # rt = int(t)
                    # if abs(int(t)-t) < 1e-5:
                    rt = int(round(t, 0)) + 1
                    # print(rt)
                    DNA[bead - self.start, rt:] += (tp1 - t) * fork_speed * dt

                # print
                # print(DNA)
            # raise
                # print(DNA[pos - self.start])

        # print(DNA[0,::])
        DNA[DNA > 1] = 1
        # print(DNA)
        # if self.number == 0:
    #        print (DNA[:10])
#            print (np.sum(DNA, axis=0)[:10])
        # raise
        self.cache = True
        self.max_t = max_t
        self.sDNA = np.sum(DNA, axis=0)
        self.DNA = DNA
        return self.sDNA, DNA

    def get_firing_time_It(self, fork_speed, normed=False, cut=0, dt=1):

        max_t = int(self.t / dt) + int(1 / fork_speed / dt) + 2

        firing_time = []
        for m in self.modules + self.ended:
            if not m.move:
                if m.activated and m.activation_time is not None:
                    firing_time.append(int(round(m.activation_time / 1.0 / dt, 0)))
        firing_time.sort()
        # print(firing_time)
        It = np.zeros(max_t)
        if cut != 0:
            for el in range(cut):
                firing_time.pop(-1)
        # print(self.t)
        for el in firing_time:
            # print(el,)
            It[int(round(el, 0))] += 1

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

    def get_free_origins_time(self, fork_speed, normed=False, dt=1):
        max_t = int(self.t / dt) + int(1 / fork_speed / dt) + 2

        free = np.zeros(max_t)

        # print("T", int(self.t), self.origins)
        for m in self.modules:
            if not m.move:
                # print(m.activation_time,m.position)

                free[:max_t - 1] += 1

        for m in self.ended:
            if not m.move:
                # print(m.activation_time,m.position)

                if m.activated:
                    free[:int(round(m.activation_time / dt, 0))] += 1
                    # print(m.activation_time,m.activated,m.passivated)
                else:
                    # Passivated origins
                    if m.activation_time is None:
                        print("Warning incorrect calculation free origins")
                    else:
                        free[:int(round(m.activation_time / dt, 0))] += 1
                        # print(m.activation_time,m.activated,m.passivated)

                    # print("Passivated")

                    # print(free)
        if normed:
            rep_p = np.array(self.get_replication_profile(fork_speed=fork_speed))
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
                m.update_position(dt=dt, fork_speed=fork_speed, mini=self.start, maxi=self.end)

        N_mod = len(self.modules)

        #####################################

        found = True
        while found:
            im = 0
            found = False
            while im < N_mod:
                # Take care of passivated Origin Left to right
                m = self.modules[im]

                if m.move:
                    if im != N_mod - 1 and m.position > self.modules[im + 1].position:
                        if self.modules[im + 1].origin:
                            passivated_origin.append(self.modules[im + 1].tag)
                            found = True
                            self.modules[im + 1].passivated = True
                            self.modules[im + 1].activation_time = self.t - \
                                dt   # To take care of 1/ dt
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

        found = True
        while found:
            im = N_mod - 1
            found = False
            while im > 0:
                # Take care of passivated Origin Right to left
                m = self.modules[im]
                if m.move:
                    if im != 0 and m.position < self.modules[im - 1].position:
                        if self.modules[im - 1].origin:
                            found = True
                            passivated_origin.append(self.modules[im - 1].tag)
                            self.modules[im - 1].passivated = True
                            self.modules[im - 1].activation_time = self.t - \
                                dt  # To take care of 1/ dt

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

        to_remove = []
        im = 0
        while im < N_mod:
            # Take care of fork collision and bond motion
            m = self.modules[im]
            if m.move:
                if im != N_mod - 1 and m.position > self.modules[im + 1].position:
                    if self.modules[im + 1].move:
                        # Collision
                        to_release.append(m.bond_tag)
                        to_release.append(self.modules[im + 1].bond_tag)
                        diff_diff.append(m.tag)
                        diff_diff.append(self.modules[im + 1].tag)
                        bind_diff.append([m.tag, self.modules[im + 1].tag])
                        outpos = (m.position + self.modules[im + 1].position) / 2.
                        if m.d == 1 and outpos < m.path[-1].pos:
                            m.path.pop(-1)
                        elif m.d == -1 and outpos > m.path[-1].pos:
                            m.path.pop(-1)
                        m.path[-1].outpos = outpos

                        if self.modules[im + 1].d == 1 and outpos < self.modules[im + 1].path[-1].pos:
                            self.modules[im + 1].path.pop(-1)
                        elif self.modules[im + 1].d == -1 and outpos > self.modules[im + 1].path[-1].pos:
                            self.modules[im + 1].path.pop(-1)

                        self.modules[im + 1].path[-1].outpos = outpos

                        to_remove.append(im)
                        to_remove.append(im + 1)
                        self.add_event(m.tag, "E")
                        self.add_event(self.modules[im + 1].tag, "E")

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
                elif m.position <= self.start or m.position >= self.end + 1:
                    # Take care of fork outside of boundaries
                    alone.append(m.tag)
                    self.add_event(m.tag, "E")
                    # print(alone,m.position)
                    """
                    if self.number == 0:
                        if  m.position > self.end:
                            print("la")
                            print(m)
                            print(self.modules)"""

                    to_release.append(m.bond_tag)
                    # print("Inside",m.position)

                    if m.position <= self.start:
                        if m.path[-1].pos <= self.start:
                            # pass
                            m.path.pop(-1)
                        else:
                            m.path[-1].outpos = self.start

                    if m.position >= self.end + 1:
                        if m.path[-1].bead >= self.end + 1:
                            # pass
                            # print(m.path, self.end+1)
                            # print("la")
                            m.path.pop(-1)
                        else:
                            m.path[-1].outpos = self.end + 1
                    """
                    if m.position > self.end + 1:
                        m.path.pop(-1)
                    elif  m.position == self.end + 1:
                        m.path[-1].outpos = self.end + 1
"""

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

            if m.path != []:
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
