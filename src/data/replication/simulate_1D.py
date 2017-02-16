import numpy as np
from replication.PMotion import Polymer


class simulate:

    def __init__(self, nori, ndiff, lengths, p_on, p_off, only_one=False, tolerance=0.1):

        self.p_on = p_on
        self.p_off = p_off
        self.ndiff = ndiff
        self.oris = []
        self.only_one = only_one
        self.lengths = lengths
        # print(nori)

        assert(len(lengths) == len(nori))

        starts = [0] + np.cumsum(lengths).tolist()
        # print(starts)

        for el, length, start in zip(nori, lengths, starts):
            if type(el) == list:
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

        # print(self.oris)
        # print(self.poly.modules)

    def simulate(self, n):

        for p in self.polys:
            p.increment_time(1)  # To start at t == 1 to check for fully replicated

        alones = 0
        for time in range(n):
            order = np.arange(self.ndiff)
            np.random.shuffle(order)
            for diff in order:
                #print(self.libre[diff],self.libre[diff] == 0)
                if self.libre[diff] == 0:
                    # print(np.random.rand())
                    if np.random.rand() < self.p_on:
                        # print("La")
                        ori_libre = []
                        for ip, p in enumerate(self.polys):
                            ori_libre += [[ip, orip] for orip in list(p.get_free_origins())]
                        if len(ori_libre) == 0:
                            continue
                        # print(ori_libre)
                        choice = np.random.randint(len(ori_libre))
                        # print("La")
                        # print(ori_libre[choice],self.poly.bound_to_origin[ori_libre[choice]])
                        # print(time)
                        # print(self.poly.state())
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

                        else:
                            if not self.only_one:
                                self.libre[diff] = 1
                                self.origins[diff] = ori_libre[choice]

                            else:
                                # If only one needed start the fork and reciclate this
                                self.polys[what_p].add_fork(
                                    [diff, diff], what_ori, [None, None], None)
                                self.libre[diff] = 2
                                self.origins[diff] = ori_libre[choice]

                elif self.libre[diff] == 1:
                    if np.random.rand() < self.p_off:
                        self.libre[diff] = 0
                        what_p, what_ori = self.origins[diff]
                        self.polys[what_p].dettach_one_diff(diff, what_ori)
                        self.origins[diff] = None

                elif self.libre[diff] == 2:
                    pass

            ended = 0
            for P in self.polys:
                bind_diff, diff_diff, update_bond, passivated_origin, to_release, alone = P.increment_time(
                    1)

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

                    for diff1 in alone:
                        self.libre[diff1] = 0
                        self.origins[diff1] = None
                else:
                    for diff1, diff2 in bind_diff:
                        for diff_t in range(self.ndiff):
                            if self.libre[diff_t] == 2:
                                self.libre[diff_t] = 0
                                self.origins[diff] = None
                                break

                    for diff in alone:
                        alones += 1

                    if alones >= 2:
                        for diff_t in range(self.ndiff):
                            if self.libre[diff_t] == 2:
                                self.libre[diff_t] = 0
                                self.origins[diff] = None
                                alones -= 2
                                break

                # Free the diff that where on origins that are now passivated
                for p in passivated_origin:
                    found = 0
                    for diff_t in range(self.ndiff):
                        if self.origins[diff_t] is not None and self.origins[diff_t][1] == p:
                            self.origins[diff_t] = None
                            self.libre[diff_t] = 0
                            found += 1
                    if found == 2:
                        print("Passivated origins with two diff")
                        raise

                """
                if self.only_one:
                    activ = 0
                    for diff_t in range(self.ndiff):
                        if self.libre[diff_t] == 2:
                            activ += 1
                    if not (activ-self.poly.get_n_activ_fork() / 2) <= 1 :
                        print(activ,self.poly.get_n_activ_fork()/2,alones)
                        raise
                """

                if P.modules == []:
                        # np.sum(np.array(P.get_replication_profile()) == 0 ) == 0:
                    ended += 1
                    # print("Ended",time)
                    #print(self.poly.get_replication_profile() == 0 )
                    # break

            if ended == len(self.lengths):
                break


if __name__ == "__main__":
    S = simu([[1, 4, 89]], 10, [1000], p_on=0.8, p_off=0.2, only_one=True)
    S.simulate(100)
    S = simu([3], 10, [1000], p_on=0.8, p_off=0.2, only_one=True)
    S.simulate(100)
