from .simulate_1D import simulate
import numpy as np
import _pickle as cPickle
from collections import namedtuple


class ensembleSim:

    def __init__(self, Nsim, Nori, Ndiff, lengths,
                 p_on, p_off, only_one, all_same_ori=False,
                 dt_speed=1,
                 fork_speed=1,
                 l_ori=[], cut=10):
        self.Nsim = Nsim
        self.Nori = Nori
        self.Ndiff = Ndiff
        self.lengths = lengths
        if type(lengths) == str:
            print("Lengths = %s" % lengths)
            raise
        self.p_on = p_on
        self.p_off = p_off
        self.only_one = only_one
        self.all_same_ori = all_same_ori
        self.dt_speed = dt_speed
        self.fork_speed = fork_speed
        self.cut = cut
        self.l_ori = l_ori

    def run_all(self, run_length=200, load_from_file=None):

        self.aIts = []
        self.aFts = []
        self.aFds = []
        self.aRps = []
        self.aDNAs = []
        self.raDNAs = []
        self.aUnrs = []
        for sim in range(self.Nsim):
            # print(sim)

            ori = self.Nori
            if self.l_ori != []:
                ori = self.l_ori

            if load_from_file is None:
                S = simulate(ori,
                             self.Ndiff,
                             self.lengths,
                             self.p_on,
                             self.p_off,
                             self.only_one,
                             dt_speed=self.dt_speed,
                             fork_speed=self.fork_speed)

                S.simulate(run_length)
            else:
                # print(sim)
                Simu = namedtuple("Simu", ["polys", "oris"])
                with open("%s%i/" % (load_from_file, sim + 1) + "polymer_timing.dat", "rb") as f:
                    polys = cPickle.load(f)
                    oris = [np.array(p.origins) - p.start for p in polys]
                    S = Simu(polys, oris)

            if sim == 0 and self.all_same_ori:
                self.l_ori = S.oris

            self.aIts.append([])
            self.aFts.append([])
            self.aFds.append([])
            self.aRps.append([])
            self.aDNAs.append([])
            self.raDNAs.append([])
            self.aUnrs.append([])

            for poly in S.polys:
                ft, it = poly.get_firing_time_It(normed=False)
                fd = poly.get_fork_density(self.cut, normed=False)
                norm = poly.get_norm()

                self.aIts[-1].append(it)
                self.aFts[-1].append(ft)
                self.aFds[-1].append(fd)
                self.aUnrs[-1].append(norm)

                self.aRps[-1].append(poly.get_replication_profile())
                self.raDNAs[-1].append(poly.get_DNA_with_time())

            unr = np.sum(np.array(self.aUnrs[-1]), axis=0)
            unr[unr == 0] = 1
            self.aIts[-1] = np.sum(np.array(self.aIts[-1]), axis=0) / unr
            self.aFds[-1] = np.sum(np.array(self.aFds[-1]), axis=0) / unr
            # print(self.raDNAs)
            self.aDNAs[-1] = 1 + np.sum(np.array(self.raDNAs[-1]), axis=0) / np.sum(self.lengths)

    def get_quant(self, name, shift=0):
        prop = getattr(self, name)
        # print(prop)
        maxl = max(map(len, prop))

        normed_prop = np.zeros((self.Nsim, maxl))
        for iIt, It in enumerate(prop):
            normed_prop[iIt, :len(It)] = It
            if shift != 0:
                normed_prop[iIt, len(It):] = It[-1]

        x = np.arange(maxl)
        y = np.mean(normed_prop, axis=0)
        err = np.std(normed_prop, axis=0)
        return x, y, err, normed_prop

    def get_rep_profile(self):
        rep = []
        for il, l in enumerate(self.lengths):
            rep.append(np.zeros(l))
            for sim in range(self.Nsim):
                rep[il] += np.array(self.aRps[sim][il]) / self.Nsim
        return rep

    def get_mean_copie(self, time):
        copie = []
        std_copie = []
        for il, l in enumerate(self.lengths):
            copie.append(np.ones((self.Nsim, l)))
            for sim in range(self.Nsim):
                copie[il][sim, np.array(self.aRps[sim][il]) < time] = 2

            std_copie.append(np.std(copie[il], axis=0))
            copie[il] = np.mean(copie[il], axis=0)

        return copie, std_copie

    def Its(self):
        return self.get_quant("aIts")

    def Fds(self):
        return self.get_quant("aFds")

    def Rps(self):
        return self.get_quant("aRps")

    def DNAs(self):
        return self.get_quant("aDNAs", shift=2)

    def n_activated_oris(self):
        return list(map(len, np.concatenate(self.aFts)))

    def error_DNA_time(self, plot=False):

        # https://academic.oup.com/nar/article/42/1/e3/2437422/The-dynamics-of-genome-replication-using-deep
        point = [(4.3714285714285808, 1.0420168067226889), (9.2571428571428562, 1.0126050420168067), (14.40000000000002, 1.0714285714285714), (17.228571428571435, 1.0420168067226889), (19.800000000000015, 0.97058823529411764), (24.428571428571431, 0.96218487394957974), (30.085714285714289, 0.97478991596638642), (32.657142857142873, 1.0714285714285714), (34.71428571428573, 1.1596638655462184), (37.028571428571425, 1.2983193277310923),
                 (39.85714285714284, 1.3277310924369747), (42.428571428571445, 1.3067226890756303), (44.48571428571428, 1.5462184873949578), (46.800000000000026, 1.588235294117647), (49.371428571428581, 1.6470588235294117), (54.771428571428551, 1.672268907563025), (59.914285714285718, 1.8613445378151261), (69.942857142857122, 1.9957983193277311), (79.971428571428589, 1.9495798319327733), (89.742857142857147, 1.8781512605042017)]
        #x_exp,y_exp = zip(*point)

        x, y, std, alls = self.DNAs()
        error = 0
        Np = 0
        shift = 0
        for xe, ye in point:
            if xe >= shift:
                i = np.argmin((x - xe + shift)**2)
                # print(x[i],xe)
                error += (ye - y[i])**2
                Np += 1
        if plot:
            return zip(*point)

        return error, Np

    def error_firing_time(self, plot=False):

        # https://academic.oup.com/nar/article/42/1/e3/2437422/The-dynamics-of-genome-replication-using-deep
        point = [(5, 0.01), (13, 0.02), (16, 0.04), (20, 0.07), (25, 0.02),
                 (30, 0.01)] + [(i, 0) for i in range(31, 70, 2)]
        unity = 1  # we want it by minutes
        point = [(time, value * unity) for time, value in point]
        x, y, std, alls = self.Its()
        error = 0
        Np = 0
        shift = 0
        for xe, ye in point:
            if xe >= shift:
                i = np.argmin((x - xe + shift)**2)
                # print(x[i],xe)
                error += (ye - y[i])**2
                Np += 1
        if plot:
            return zip(*point)
        return error, Np
