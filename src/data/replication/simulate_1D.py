import numpy as np
import copy
from replication.PMotion import Polymer, Diffusing


class simulate:

    def __init__(self, nori, ndiff, lengths, p_on, p_off, only_one=False,
                 fork_speed=1, dt_speed=1, tolerance=0.1,
                 gindin=True, p_v=1, random=False, positions=None, ramp=None, max_ramp=None,
                 ramp_type="linear", strengths=None):

        self.p_on = p_on
        self.p_off = p_off
        self.ndiff = ndiff
        self.oris = []
        self.only_one = only_one
        self.lengths = lengths
        self.dt_speed = dt_speed
        self.fork_speed = fork_speed
        self.gindin = gindin
        self.p_v = p_v
        self.Ndiff_libre_t = []
        self.positions = positions
        self.ramp = ramp
        self.max_ramp = max_ramp
        self.strengths = copy.deepcopy(strengths)
        self.uniform_strength = False
        self.ramp_type = ramp_type
        if self.ramp_type not in ["linear", "exp", "pulse"]:
            raise "Invalid ramp type"

        if self.max_ramp and self.max_ramp > self.ndiff:
            print(self.max_ramp, self.ndiff)
            raise "Problem max_ramp > ndiff"

        # print(nori)

        assert(len(lengths) == len(nori))

        starts = [0] + np.cumsum(lengths).tolist()
        # print(starts)

        for el, length, start in zip(nori, lengths, starts):
            if type(el) in [list, np.ndarray]:
                self.oris.append(el)
                continue
            self.oris.append([])

            while (el - len(self.oris[-1])) / el > tolerance:
                self.oris[-1] = [np.random.randint(length) for ori in range(el)]
                self.oris[-1] = list(set(self.oris[-1]))
        # print(self.oris)
            # print(len(self.oris),(nori-len(self.oris)) / nori)
        if strengths is None or strengths == []:
            self.uniform_strength = True
            self.strengths = [[1] * len(o) for o in self.oris]

        if positions is None:
            self.polys = [Polymer(i, start, start + length - 1, np.array(oris) + start, random=random, strengths=strengths) for i, (start, length, oris, strengths) in
                          enumerate(zip(starts, lengths, self.oris, self.strengths))]
        else:
            self.polys = [Polymer(i, start, start + length - 1, np.array(oris) + start, random=random, positions=np.array(position) + start, strengths=strengths) for i, (start, length, oris, position, strengths) in
                          enumerate(zip(starts, lengths, self.oris, self.positions, self.strengths))]

        class MyD(dict):

            def __init__(self, *args):
                dict.__init__(self, *args)
                self.set = 0

            def __setitem__(self, key, value):
                if key == 0 and value == 0:
                    self.set += 1
                if self.set == 1:
                    raise
                dict.__setitem__(self, key, value)

        self.libre = {d: 0 for d in np.arange(ndiff)}
        self.origins = {d: None for d in np.arange(ndiff)}
        self.record_diffusing = [Diffusing(d) for d in np.arange(ndiff)]

        # print(self.oris)
        # print(self.poly.modules)

    def get_free(self, maxt=None):
        V = []
        maxt = self.time + int(1 / self.fork_speed) + 2
        # print(maxt)
        for k in self.record_diffusing:
            # print(k)
            V.append(np.array(k.build_time_line(maxt=maxt)))
        V = np.array(V)
        return np.sum((V == 0) + (V == 1), axis=0)

    def simulate(self, n):
        """
        for P in self.polys:
            # events =
            bind_diff, diff_diff, update_bond, passivated_origin, to_release, alone = P.increment_time(
                dt=self.dt_speed / 2, fork_speed=self.fork_speed)
        """
        alones = []
        extra = 0

        ori_libre = []
        dt_speed_d2 = self.dt_speed  # / 2.

        for ip, p in enumerate(self.polys):
            ori_libre += [[ip, orip] for orip in list(p.get_free_origins())]

        for time in range(n):

            #print(time // 2, len(ori_libre))
            ended = 0
            #print("Start", len(ori_libre))
            if time != 0:
                #print("Mv fork")
                for P in self.polys:
                    # events =
                    bind_diff, diff_diff, update_bond, passivated_origin, to_release, alone = P.increment_time(
                        dt=dt_speed_d2)
                    if self.only_one:
                        for k in P.bound_to_origin:
                            if P.bound_to_origin[k] != []:
                                print("found an alone diffusing element that should have started")
                                raise

                    if not self.only_one:
                        for diff1, diff2 in bind_diff:
                            self.libre[diff1] = 0
                            self.libre[diff2] = 0
                            self.origins[diff1] = None
                            self.origins[diff2] = None
                            self.record_diffusing[diff1].end_replication(time * dt_speed_d2)
                            self.record_diffusing[diff2].end_replication(time * dt_speed_d2)

                        for diff1 in alone:
                            # or it was at the end of a chromosome or on an origin passivated
                            if self.libre[diff1] == 2:
                                self.record_diffusing[diff1].end_replication(time * dt_speed_d2)
                            elif self.libre[diff1] == 1:
                                self.record_diffusing[diff1].end_bound(time * dt_speed_d2)
                            else:
                                print("Very strange (simulate_1D)")
                                raise

                            self.libre[diff1] = 0
                            self.origins[diff1] = None
                            """
                            if self.record_diffusing[diff1].replicating != []:
                                if len(self.record_diffusing[diff1].replicating[-1]) == 2:
                                    pass
                                else:
                                    self.record_diffusing[diff1].end_replication(
                                        time * dt_speed_d2)
                            """
                    else:

                        for diff1, diff2 in bind_diff:
                            # Release one randoms occupied diffusing element
                            diff = diff1 if diff1 is not None else diff2
                            if diff1 is not None and diff2 is not None:
                                print("Strange")
                                raise

                            if diff is None:
                                for diff_t in range(self.ndiff):
                                    if self.libre[diff_t] == 2:
                                        self.libre[diff_t] = 0
                                        self.origins[diff_t] = None
                                        break
                            else:
                                if self.libre[diff] == 2:
                                    self.libre[diff] = 0
                                    self.origins[diff] = None
                                    self.record_diffusing[diff].end_replication(
                                        time * dt_speed_d2)

                                else:
                                    print("Problem")
                                    raise

                        for diff in alone:
                            alones.append(diff)

                        if len(alones) >= 2:
                            isnone = np.equal(alones, None)
                            if len(isnone) == 0 or len(isnone) == len(alones):
                                pass
                            else:

                                diff = np.array(alones)[~isnone][0]

                                self.libre[diff] = 0
                                self.origins[diff] = None
                                self.record_diffusing[diff].end_replication(time * dt_speed_d2)
                                diff1 = alones.pop(alones.index(diff))
                                diff2 = alones.pop(alones.index(None))
                                print("Freed", diff1, diff2)
                                if diff1 is not None and diff2 is not None:
                                    print(diff1, diff2, alones)
                                    print("Strange")
                                    raise

                    # Free the diff that where on origins that are now passivated

                    for p in passivated_origin:
                        found = 0
                        for diff_t in range(self.ndiff):
                            if self.origins[diff_t] is not None and self.origins[diff_t][1] == p:
                                self.origins[diff_t] = None
                                self.libre[diff_t] = 0
                                self.record_diffusing[diff_t].end_bound(time * dt_speed_d2)
                                found += 1
                        if found == 2:
                            print("Passivated origins with two diff")
                            raise

                        for iori, (ip, ori) in enumerate(ori_libre):
                            found = False
                            if ori == p:
                                ori_libre.pop(iori)
                                found = True
                                break
                        if not found:
                            print(p)
                            print(ori_libre)
                            print("Missing origin")
                            raise
                    if P.modules == []:
                            # np.sum(np.array(P.get_replication_profile()) == 0 ) == 0:
                        ended += 1

            # Random checks
            if np.random.random() < 0.01:
                ori1_libre = []
                for ip, p in enumerate(self.polys):
                    ori1_libre += [[ip, orip] for orip in list(p.get_free_origins())]
                try:
                    assert(len(ori1_libre) == len(ori_libre))
                except:
                    print(ori_libre)
                    print(ori1_libre)
                    print(self.polys[0].ended)

                    raise
            # print(time)
            # if time % 2 == 1:
            #    continue
            # print("starting")
            order = np.arange(self.ndiff)

            if self.ramp:
                order = int(self.ramp * time * self.dt_speed)
                if self.max_ramp:
                    if self.ramp_type == "linear":
                        order = min(order, self.max_ramp)
                    if self.ramp_type == "exp":
                        order = self.max_ramp * (1 - np.exp(- time * self.dt_speed / self.ramp))

                    if self.ramp_type == "pulse":
                        order = self.max_ramp
                order = np.arange(int(order))

            np.random.shuffle(order)
            self.Ndiff_libre_t.append(
                np.sum([1 for tmp_diff in order if self.libre[tmp_diff] == 0]))
            #Nori_libre = len(ori_libre)
            # print("Middle")
            # print(len(ori_libre))
            #print("Start ori", len(ori_libre))
            for diff in order:
                # print(self.libre[diff],self.libre[diff] == 0)
                if self.libre[diff] == 0:
                    if len(ori_libre) == 0:
                        continue

                    # print(np.random.rand())
                    if self.gindin:
                        if np.random.rand() < self.p_on:
                            # Interaction

                            pass

                        else:
                            continue

                    elif not self.gindin:
                        Nori_libre = len(ori_libre)

                        """
                        n_possible_ori = np.random.binomial(Nori_libre, self.p_v)
                        if n_possible_ori >= 1:
                            self.record_diffusing[diff].in_volume_of_interaction(
                                [ori for ori in range(n_possible_ori)], time * dt_speed_d2)

                        p_one_interaction = 1 - (1 - self.p_on)**n_possible_ori
                        if np.random.rand() < p_one_interaction:
                            # Interaction
                            pass
                        else:
                            continue
                        """
                        if self.uniform_strength:
                            ln_no_interaction = Nori_libre * np.log(1 - self.p_on * self.p_v)
                            if np.log(np.random.rand()) > ln_no_interaction:
                                # Interaction
                                pass
                            else:
                                continue
                        else:
                            ln_no_interaction = Nori_libre * np.log(1 - self.p_v)
                            if np.log(np.random.rand()) > ln_no_interaction:
                                # Interaction
                                pass
                            else:
                                continue

                    choice = np.random.randint(len(ori_libre))
                    what_p, what_ori = ori_libre[choice]

                    if not self.uniform_strength:
                        if not (np.random.rand() < self.p_on * self.polys[what_p].o_strength[what_ori]):
                            continue

                    two = self.polys[what_p].attach_one_diff(diff, what_ori, None)
                    if two:
                        # print("Start",ori_libre[choice])
                        # print(self.poly.modules)
                        [diff1, _], [diff2, _] = self.polys[what_p].get_diff_at_origin(what_ori)
                        self.polys[what_p].add_fork(
                            [diff1, diff2], what_ori, [None, None], None, fork_speed=self.fork_speed)
                        self.libre[diff1] = 2
                        self.libre[diff2] = 2
                        if diff == diff1:
                            self.record_diffusing[diff2].end_bound(time * dt_speed_d2)
                        else:

                            self.record_diffusing[diff1].end_bound([time * dt_speed_d2, "TH"])

                        self.record_diffusing[diff1].start_replication(
                            what_ori, time * dt_speed_d2)
                        self.record_diffusing[diff2].start_replication(
                            what_ori, time * dt_speed_d2)

                        ori_libre.remove(ori_libre[choice])

                    else:
                        if not self.only_one:
                            self.libre[diff] = 1
                            self.origins[diff] = ori_libre[choice]
                            self.record_diffusing[diff].start_bound(what_ori, time * dt_speed_d2)

                        else:
                            # If only one needed start the fork and reciclate this
                            self.polys[what_p].add_fork(
                                [diff, None], what_ori, [None, None], None, fork_speed=self.fork_speed)
                            self.libre[diff] = 2
                            self.origins[diff] = ori_libre[choice]
                            self.record_diffusing[diff].start_replication(
                                what_ori, time * dt_speed_d2)

                            ori_libre.remove(ori_libre[choice])

                elif self.libre[diff] == 1:
                    if np.random.rand() < self.p_off:
                        self.libre[diff] = 0
                        what_p, what_ori = self.origins[diff]
                        self.polys[what_p].dettach_one_diff(diff, what_ori)
                        self.origins[diff] = None
                        self.record_diffusing[diff].end_bound(time * dt_speed_d2)

                elif self.libre[diff] == 2:
                    pass

            if ended == len(self.lengths):
                extra += 1
                self.time = time * dt_speed_d2

                if extra > 10:
                    # print(extra)
                    break


if __name__ == "__main__":
    S = simulate([[1, 4, 89]], 10, [1000], p_on=0.8, p_off=0.2, only_one=True)
    S.simulate(100)
    S = simulate([3], 10, [1000], p_on=0.8, p_off=0.2, only_one=True)
    S.simulate(100)
