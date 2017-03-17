import numpy as np
from replication.PMotion import Polymer, Diffusing


class simulate:

    def __init__(self, nori, ndiff, lengths, p_on, p_off, only_one=False,
                 fork_speed=1, dt_speed=1, tolerance=0.1, gindin=True, p_v=1):

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
            #print(len(self.oris),(nori-len(self.oris)) / nori)
        self.polys = [Polymer(i, start, start + length - 1, np.array(oris) + start) for i, (start, length, oris) in
                      enumerate(zip(starts, lengths, self.oris))]

        self.libre = {d: 0 for d in np.arange(ndiff)}
        self.origins = {d: None for d in np.arange(ndiff)}
        self.record_diffusing = [Diffusing(d) for d in np.arange(ndiff)]


        # print(self.oris)
        # print(self.poly.modules)

    def get_free(self, maxt=None):
        V = []
        for k in S.record_diffusing:
            V.append(np.array(k.build_time_line(maxt=maxt)))
        return np.sum( (V==0) + ( V == 1),axis=0)

    def simulate(self, n):

        alones = []

        ori_libre = []
        for ip, p in enumerate(self.polys):
            ori_libre += [[ip, orip] for orip in list(p.get_free_origins())]

        for time in range(n):

            ended = 0
            if time != 0:
                for P in self.polys:
                    bind_diff, diff_diff, update_bond, passivated_origin, to_release, alone = P.increment_time(
                        dt=self.dt_speed, fork_speed=self.fork_speed)
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
                            self.record_diffusing[diff1].end_replication(time * self.dt_speed)
                            self.record_diffusing[diff2].end_replication(time * self.dt_speed)


                        for diff1 in alone:
                            # or it was at the end of a chromosome or on an origin passivated
                            self.libre[diff1] = 0
                            self.origins[diff1] = None
                            if self.record_diffusing[diff1].replicating != []:
                                if len(self.record_diffusing[diff1].replicating[-1]) == 2:
                                    pass
                                else:
                                    self.record_diffusing[diff1].end_replication(time * self.dt_speed)

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
                                    self.record_diffusing[diff].end_replication(time * self.dt_speed)

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
                                self.record_diffusing[diff].end_replication(time * self.dt_speed)
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
                                self.record_diffusing[diff_t].end_bound(time * self.dt_speed)
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
                        # print("Ended",time)
                        #print(self.poly.get_replication_profile() == 0 )
                        # break

            if ended == len(self.lengths):
                self.time = time
                break

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


            order = np.arange(self.ndiff)
            np.random.shuffle(order)

            for diff in order:
                #print(self.libre[diff],self.libre[diff] == 0)
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
                        n_possible_ori = np.random.binomial(Nori_libre, self.p_v)
                        if n_possible_ori >= 1:
                            self.record_diffusing[diff].in_volume_of_interaction([ori for ori in range(n_possible_ori)], time * self.dt_speed)

                        p_one_interaction = 1 - (1 - self.p_on)**n_possible_ori
                        if np.random.rand() < p_one_interaction:
                            # Interaction
                            pass
                        else:
                            continue

                    choice = np.random.randint(len(ori_libre))
                    what_p, what_ori = ori_libre[choice]

                    two = self.polys[what_p].attach_one_diff(diff, what_ori, None)
                    if two:
                        # print("Start",ori_libre[choice])
                        # print(self.poly.modules)
                        [diff1, _], [diff2, _] = self.polys[what_p].get_diff_at_origin(what_ori)
                        self.polys[what_p].add_fork(
                            [diff1, diff2], what_ori, [None, None], None)
                        self.libre[diff1] = 2
                        self.libre[diff2] = 2
                        self.record_diffusing[diff1].start_replication(what_ori, time * self.dt_speed)
                        self.record_diffusing[diff2].start_replication(what_ori, time * self.dt_speed)

                        ori_libre.remove(ori_libre[choice])

                    else:
                        if not self.only_one:
                            self.libre[diff] = 1
                            self.origins[diff] = ori_libre[choice]
                            self.record_diffusing[diff].start_bound(what_ori, time * self.dt_speed)

                        else:
                            # If only one needed start the fork and reciclate this
                            self.polys[what_p].add_fork(
                                [diff, None], what_ori, [None, None], None)
                            self.libre[diff] = 2
                            self.origins[diff] = ori_libre[choice]
                            self.record_diffusing[diff].start_replication(what_ori, time * self.dt_speed)

                            ori_libre.remove(ori_libre[choice])

                elif self.libre[diff] == 1:
                    if np.random.rand() < self.p_off:
                        self.libre[diff] = 0
                        what_p, what_ori = self.origins[diff]
                        self.polys[what_p].dettach_one_diff(diff, what_ori)
                        self.origins[diff] = None
                        self.record_diffusing[diff].end_bound(time * self.dt_speed)


                elif self.libre[diff] == 2:
                    pass


if __name__ == "__main__":
    S = simulate([[1, 4, 89]], 10, [1000], p_on=0.8, p_off=0.2, only_one=True)
    S.simulate(100)
    S = simulate([3], 10, [1000], p_on=0.8, p_off=0.2, only_one=True)
    S.simulate(100)
