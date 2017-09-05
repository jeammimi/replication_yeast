from .simulate_1D import simulate
import numpy as np
import _pickle as cPickle
from collections import namedtuple
import os
from tqdm import tqdm
import pandas as pd
import h5py
import json
from scipy.stats import poisson
from replication.tools import load_ori_position, load_lengths_and_centro


class ensembleSim:

    def __init__(self, Nsim, Nori, Ndiff, lengths,
                 p_on, p_off, only_one, all_same_ori=True,
                 dt_speed=1,
                 fork_speed=1,
                 gindin=True,
                 p_v=1,
                 l_ori=[], cut=10, random=False, one_minute=False,
                 positions=None, ramp=None,
                 max_ramp=None, ramp_type="linear", strengths=[], hdf5_file=None,
                 D_Ndiff="pulse"):
        self.Nsim = Nsim
        self.Nori = Nori
        self.Ndiff = Ndiff
        self.lengths = lengths
        if type(lengths) == str:
            print("Lengths = %s" % lengths)
            raise
        if lengths and type(lengths[0]) == list:
            print("lengts = ", lengths)
            print("But should be a list")
            raise
        assert(type(gindin) == bool)
        assert(type(only_one) == bool)
        self.p_on = p_on
        self.p_off = p_off
        self.only_one = only_one
        self.all_same_ori = all_same_ori
        self.dt_speed = dt_speed
        self.fork_speed = fork_speed
        self.gindin = gindin
        self.p_v = p_v
        self.cut = cut
        self.l_ori = l_ori
        self.random = random
        self.one_minute = one_minute
        self.positions = positions
        self.ramp = ramp
        self.max_ramp = max_ramp
        self.ramp_type = ramp_type
        self.strengths = strengths
        self.hdf5_file = None
        self.D_Ndiff = D_Ndiff

    def add_precomputed(self, name, file_hdf5="None", precision=None, two=False):
        qt = getattr(self, name)()
        with h5py.File(file_hdf5, 'a') as myfile:

            quant = myfile.get("analysis")
            if myfile.get("analysis") is None:
                quant = myfile.create_group("analysis")

            if quant.get(name) is not None:
                print(name, "Allready computed")
                return
            # print(quant.get(name))
            # print(type(qt[0]))
            if qt != [] and type(qt[0]) in[list, np.ndarray]:

                prop = quant.create_group(name)
                if precision:
                    prop.create_dataset("precision", data=precision)
                maxi = None
                if two:
                    maxi = 2
                for i in range(len(qt[:maxi])):
                    if precision:
                        prop.create_dataset(str(i), data=list(
                            map(lambda x: int(x * precision), qt[i])))
                    else:
                        prop.create_dataset(str(i), data=np.array(qt[i]))
            else:
                prop = quant.create_dataset(name, data=qt)

    def show_parameters(self, show_ori=True):
        P = ["Nsim", "Nori", "Ndiff", "lengths", "p_on", "p_off",
             "only_one", "all_same_ori", "dt_speed",
             "fork_speed", "gindin", "p_v", "cut", "l_ori", "ramp", "max_ramp"]

        for parameter in P:
            if (parameter == "l_ori" or parameter == "Nori") and not show_ori:
                print(parameter, self.nori)
                continue
            if hasattr(self, parameter):
                print(parameter, getattr(self, parameter))
            else:
                print(parameter, "Not defined")

    def data(self):

        return [self.aIts,
                self.aFts,
                self.aFds,
                self.aRps,
                self.aDNAs,
                self.raDNAs,
                self.aUnrs,
                self.aFree_origins]

    def load_data(self, data):
        self.aIts, self.aFts, self.aFds, self.aRps, self.aDNAs, self.raDNAs, self.aUnrs, self.aFree_origins = data
        unr = np.sum(np.array(self.aUnrs), axis=1)
        self.anIts = self.aIts * unr

    def remove_correlations(self):
        del self.aIODs
        del self.aIRTDs
        del self.aTLs

    def add_traj(self, N, run_length=10000):
        old_nsim = 0 + self.Nsim
        self.Nsim = N
        self.run_all(init=False)
        self.Nsim = old_nsim + N

    def run_all(self, run_length=200, load_from_file=None, correlation=True, skip=[], single=False, init=True):

        if init:
            self.aIts = []
            self.aIfs = []

            self.aFts = []
            self.aFds = []
            self.aRps = []
            self.aDNAs = []
            self.raDNAs = []
            self.aUnrs = []
            self.aFree_origins = []
            self.aFree_Diff_bis = []
            self.anIts = []
            self.aFree_Diff = []
            self.aFiring_Position = []
            self.aIODs = []
            self.aIRTDs = []
            self.aTLs = []
            self.record_diffusing = []

        found = 0
        for sim in tqdm(range(self.Nsim)):

            ori = self.Nori
            if self.l_ori != []:
                ori = self.l_ori

            # check dimension of position
            positions = self.positions
            if self.positions and type(self.positions[0][0]) is list:
                positions = self.positions[sim]

            Nd = self.Ndiff
            max_ramp = self.max_ramp
            if self.D_Ndiff == "poisson":
                Nd = poisson.rvs(size=1, mu=self.Ndiff)[0]
                max_ramp = Nd

            if load_from_file is None:
                S = simulate(ori,
                             Nd,
                             self.lengths,
                             self.p_on,
                             self.p_off,
                             self.only_one,
                             dt_speed=self.dt_speed,
                             fork_speed=self.fork_speed,
                             gindin=self.gindin,
                             p_v=self.p_v,
                             random=self.random,
                             positions=positions,
                             ramp=self.ramp,
                             max_ramp=max_ramp,
                             ramp_type=self.ramp_type,
                             strengths=self.strengths
                             )

                S.simulate(run_length)
                found += 1
                self.record_diffusing.append(S.record_diffusing)
            else:
                # print("Sim", sim)
                if sim in skip:
                    # print("skip", skip)
                    continue
                # print(sim)
                Simu = namedtuple("Simu", ["polys", "oris", "Ndiff_libre_t", "record_diffusing"])
                troot = "%s%i/" % (load_from_file, sim + 1)
                if single:
                    troot = load_from_file
                file_to_open = troot + "polymer_timing.dat"
                try:
                    if os.path.exists(file_to_open):
                        with open(file_to_open, "rb") as f:
                            polys = cPickle.load(f)
                            oris = [np.array(p.origins) - p.start for p in polys]
                        Ndiff_libre_t = []
                        if os.path.exists(troot + "Ndiff_libre_t.dat"):
                            with open(troot + "Ndiff_libre_t.dat", "rb") as f:
                                Ndiff_libre_t = cPickle.load(f)
                        record_diffusing = []
                        if os.path.exists(troot + "record_diffusing.dat"):
                            with open(troot + "record_diffusing.dat", "rb") as f:
                                record_diffusing = cPickle.load(f)
                        self.record_diffusing.append(record_diffusing)

                        S = Simu(polys, oris, Ndiff_libre_t, record_diffusing)
                        found += 1
                    else:
                        print(file_to_open, "does not exist")
                        continue
                except EOFError:
                    print("Not all files in %i readable" % sim)

            if found == 1 and self.all_same_ori:
                self.l_ori = S.oris

            self.aIts.append([])
            self.aIfs.append([])
            self.anIts.append([])
            self.aFts.append([])
            self.aFds.append([])
            self.aRps.append([])
            self.aDNAs.append([])
            self.raDNAs.append([])
            self.aUnrs.append([])
            self.aFree_Diff.append([])
            self.aFree_origins.append([])
            self.aFree_Diff_bis.append([])
            self.aFiring_Position.append([])
            self.aIODs.append([])
            self.aIRTDs.append([])
            self.aTLs.append([])

            for poly in S.polys:
                if self.one_minute:
                    dt = 1
                else:
                    dt = self.dt_speed
                # Cut == 0 because we removed them from all the chromosomes
                ft, it = poly.get_firing_time_It(
                    fork_speed=self.fork_speed, cut=0, normed=False, dt=dt)
                fd = poly.get_fork_density(fork_speed=self.fork_speed,
                                           cut=0, normed=False, dt=dt)  # Normed afteward

                self.aIts[-1].append(it)
                self.aFts[-1].append(ft)
                self.aFds[-1].append(fd)

                self.aRps[-1].append(poly.get_replication_profile(fork_speed=self.fork_speed, dt=dt))
                self.raDNAs[-1].append(poly.get_DNA_with_time(fork_speed=self.fork_speed, dt=dt)[0])
                if correlation:
                    iods, irtds, tls = poly.get_correlations(
                        fork_speed=self.fork_speed, dt=dt, thresh=0.99)
                    self.aIODs[-1].append(iods)
                    self.aIRTDs[-1].append(irtds)
                    self.aTLs[-1].append(tls)

                """
                All the following line to be able to compute No(t-1)
                """
                # print(self.aUnrs[-1][-1])
                # .append(poly.get_DNA_with_time(fork_speed=self.fork_speed)[0])
                # print(self.raDNAs[-1][-1][-1])

                Free_o = poly.get_free_origins_time(
                    fork_speed=self.fork_speed, normed=False, dt=dt).tolist()
                assert (Free_o[-1] == 0)

                self.aFree_origins[-1].append(np.array([len(poly.origins)] + Free_o[:-1]))
                # self.aFree_origins[-1].append(Free_o)
                # print(self.aFree_origins[-1])
                # assert(1 == 0)

                """
                len_poly = poly.end + 1 - poly.start

                assert(self.raDNAs[-1][-1][-1] == len_poly)

                self.raDNAs[-1][-1] = self.raDNAs[-1][-1].tolist()
                self.raDNAs[-1][-1].pop(0)
                self.raDNAs[-1][-1].append(len_poly)
                self.raDNAs[-1][-1] = np.array(self.raDNAs[-1][-1])

                # print(self.raDNAs[-1][-1])
                # self.aUnrs[-1][-1] = self.aUnrs[-1][-1]
"""
                len_poly = poly.end + 1 - poly.start

                self.aUnrs[-1].append(len_poly - self.raDNAs[-1][-1])

                ftime, firing_position = poly.get_dist_between_activated_origins(
                    fork_speed=self.fork_speed, dt=dt)
                self.aFiring_Position[-1].append(firing_position)

                # print (norm.shape,self.aUnrs[-1][-1].shape)

                # raise

                # print(it)
            DNA_time = np.sum(np.array(self.raDNAs[-1]), axis=0) / np.sum(self.lengths)
            try:
                for t in range(len(DNA_time)):
                    tp = int(round(t * dt / self.dt_speed, 0))
                    if tp > len(S.Ndiff_libre_t) - 1:
                        break
                    self.aFree_Diff_bis[-1].append(S.Ndiff_libre_t[tp])
            except:
                # Not available in 3D
                pass
            """
            try:
                self.aFree_Diff[-1] = S.get_free()
                # print(self.aFree_Diff[-1])
            except:
                pass"""

            bins = 100
            for poly in S.polys:
                self.aIfs[-1].append(poly.get_firing_at_fraction(DNA_time=DNA_time,
                                                                 fork_speed=self.fork_speed, cut=0, bins=bins))

            self.aIfs[-1] = np.sum(np.array(self.aIfs[-1]), axis=0) / \
                (np.array(np.arange(0, 1, 1 / bins) + 1 / 100.) * self.length)[::-1]
            # print (np.array(np.arange(0,1,1/bins) * np.sum(self.lengths))[::-1])

            unr = np.sum(np.array(self.aUnrs[-1]), axis=0)
            unr[unr == 0] = np.nan
            self.anIts[-1] = np.sum(np.array(self.aIts[-1]), axis=0)

            self.aIts[-1] = np.sum(np.array(self.aIts[-1]), axis=0) / unr
            self.aFds[-1] = np.sum(np.array(self.aFds[-1]), axis=0) / self.length
            self.aFree_origins[-1] = np.sum(np.array(self.aFree_origins[-1]), axis=0)
            # print(self.raDNAs)
            self.aDNAs[-1] = 1 + np.sum(np.array(self.raDNAs[-1]), axis=0) / self.length
        return S

    def get_what(self, what, fraction=[0, 1], max_track_length=None):
        """return an array which contain a concatenation by sim
        for each sim it is an array which contain a list of the given quantity for evey time step
        IOD, IRTD, or TL
        """

        def recompute(what, tl, max_track_length):
            res = []
            for ich, ch in enumerate(what):
                res.append([])

                for ipos, spos in enumerate(ch):
                    # Go throug time
                    # print(spos)
                    # print(spos,)

                    if type(spos) is not list:
                        spos = [] + spos.tolist()
                    else:
                        spos = [] + spos
                    if spos == []:
                        res[-1].append([])
                        continue

                    spos.insert(0, 0)
                    pos = np.cumsum(spos)
                    # print(tl[ich][ipos])
                    keep = np.array(tl[ich][ipos]) < max_track_length
                    kpos = pos[np.array(keep, np.bool)]
                    pos = kpos[1:] - kpos[:-1]
                    res[-1].append(pos)

                    """
                    if np.any(keep == False):
                        print(pos.shape, keep.shape, pos[keep].shape)
                        print(len(res[-1][-1]), len(ch[ipos]))
                        # print(spos, pos, keep, tl[ich][ipos])
                        print(res[-1][-1])

                        raise"""
                    # return
            return np.array(res).T

        iod3 = []
        for sim in range(self.Nsim):
            def get_by_time(what=what):
                # print(sim)
                iods = np.array(getattr(self, "a" + what + "s")[sim])
                if max_track_length is not None:
                    tl = np.array(getattr(self, "aTLs")[sim])
                    tl = tl.T
                iods = iods.T
                iods2 = []
                fraction_time = np.array(self.raDNAs[sim]).copy()

                for ichl, chl in enumerate(self.lengths):
                    # Normalise to 1 by dividing by chromosome length
                    fraction_time[ichl] /= chl

                to_keep = iods

                if max_track_length is not None:
                    # print(tl[ich].shape)
                    to_keep = recompute(iods.T, tl.T, max_track_length)

                # print(fraction_time.shape)
                for ich, (ch_what, ch_fraction) in enumerate(zip(to_keep, fraction_time.T)):
                    # We go throug time and
                    # By chromosomes select where they match the selected fraction:
                    select = (ch_fraction >= fraction[0]) * (ch_fraction <= fraction[1])
                    # print(select)
                    # return
                    if np.sum(select) >= 2:

                        iods2.append(np.concatenate(ch_what[select]))

                    if np.sum(select) == 1:
                        # print(ch_what)
                        iods2.append(np.array(ch_what[select][0]))
                        """
                        print(iods2[-1])
                        print(iods2[-2])
                        print(np.concatenate([[], []]).shape)
                        print(np.array([]).shape)
                        return"""
                    if np.sum(select) == 0:
                        iods2.append(np.array([]))
                return iods2
            iod3 += get_by_time()
        return iod3

    def get_cum_sum_hist(self, what, bins=100, fraction=[0, 1], max_track_length=None):
        """Cumulative histogram in a combing like fashion
        as the time steps are all used and added together"""

        if what != "ori":
            data = self.get_what(what, fraction=fraction, max_track_length=max_track_length)
        elif what == "ori":
            data = [np.array(io)[1:] - np.array(io)[:-1] for io in self.l_ori]

        m = []
        for i in data:
            m += i.tolist()  # np.mean(i) for i in iod3 if i != [] ]
        self.m = m
        y, x = np.histogram(m, bins=bins, normed=True)
        # hist(m,bins=100,normed=True,cumulative=-1,histtype='step')

        y = np.array([0] + np.cumsum(y).tolist())

        y /= y[-1]
        # print(y[0], y[-1])

        y = 1 - y
        # plot( 5*(x[1:]/2+x[:-1]/2),y)
        return x, y

    def get_quant(self, name, shift=0, n_rep=None, cut=0):
        if shift != 0:
            print("You should not use it")
        prop = getattr(self, name)
        # print(prop)

        times = self.get_times_replication(n_rep=n_rep)
        # print(times)
        # print(maxl)
        if -1 in times:
            maxl = int(max(map(len, prop)))
        else:
            maxl = int(max(times / self.dte))
        if name == "aIfs":
            maxl = len(prop[0])

        normed_prop = np.zeros((len(prop[:n_rep]), maxl))
        # print("Nan")
        normed_prop += np.nan

        for iIt, It in enumerate(prop[:n_rep]):
            # print(len(It), maxl)
            normed_prop[iIt, :min(len(It), maxl)] = np.array(It[:min(len(It), maxl)])

            if cut != 0 and name in ["anIts", "aFds"]:
                # Remove last cut:
                # print("Before", normed_prop[iIt])
                # print("la")
                removed = 0
                if cut != 0:
                    for i in range(1, len(normed_prop[iIt])):
                        while removed != cut and normed_prop[iIt][-i] > 0:
                            # print(i)
                            normed_prop[iIt][-i] = -1
                            removed += 1

                        if removed == cut:
                            normed_prop[iIt][-i:] = np.nan
                            break

            # print("After", normed_prop[iIt])

            if shift != 0:
                normed_prop[iIt, len(It):] = It[-1]
        self.all = normed_prop
        x = np.arange(maxl)
        if n_rep:
            y = np.nanmean(normed_prop[:n_rep], axis=0)
            err = np.std(normed_prop[:n_rep], axis=0)
        else:
            y = np.nanmean(normed_prop, axis=0)
            err = np.std(normed_prop, axis=0)

        return x * self.dte, y, err, normed_prop

    def get_times_replication(self, finished=True, n_rep=None):
        v = self.try_load_property("get_times_replication")
        if v is not None:
            return v

        times = []
        for rep in self.aRps[:n_rep]:
            times.append(-1)
            for c in rep:
                if finished and np.sum(np.equal(c, None)) != 0:
                    times[-1] = -1
                    break
                else:
                    times[-1] = max(times[-1], max(np.array(c)[~np.equal(c, None)]))

        return np.array(times) * self.dte

    @property
    def nori(self):
        nori = 1.0 * np.sum(list(map(len, self.l_ori)))
        if nori == 0:
            print("Warning, no origins ")
        return nori

    @property
    def length(self):
        return np.sum(self.lengths)

    @property
    def dte(self):
        if self.one_minute:
            return 1
        else:
            return self.dt_speed

    def try_load_property(self, name):
        # print(name)
        if hasattr(self, "hdf5_file") and self.hdf5_file is not None:
            with h5py.File(self.hdf5_file, 'r') as myfile:

                quant = myfile.get("analysis")
                if quant is not None:
                    prop = quant.get(name)
                    # print(prop, hasattr(prop, "shape"))

                    if hasattr(prop, "shape"):
                        return prop.value
                    # print(prop, dir(prop))
                    if prop is not None:
                        return [prop[str(i)].value for i in range(len(prop))]
        return None

    def get_dist_between_activated_origins(self, time=None):
        """Time in minutes"""

        v = self.try_load_property("get_dist_between_activated_origins")
        if v is not None:
            return v

        Dist = []
        if time is None:
            time = 1e8
        else:
            time = time / self.dte
        # print(time)
        for fps in self.aFiring_Position:
            for fp in fps:

                fired = fp[::, 0] <= time
                dist = fp[fired][::, 1]
                dist = dist[1:] - dist[:-1]
                Dist.extend(dist)
        return Dist

    def get_time_at_fraction(self, frac=1, bead=True):
        dna = frac + 1
        x, DNA = self.DNAs()[:2]
        # print(DNA)
        for iid, d in enumerate(DNA):
            if d >= dna:
                return x[iid]
        return x[-1]

    def It_Mean_field_origins(self, n_rep=None):
        v = self.try_load_property("It_Mean_field_origins")
        if v is not None:
            return v

        x, y = self.Free_Diff_bis(n_rep=n_rep)[:2]

        x, y1 = self.Free_origins(n_rep=n_rep)[:2]

        x, DNA = self.DNAs(n_rep=n_rep)[:2]

        Unr = (2 - DNA) * self.length

        return x, y * y1 / Unr * self.p_on * self.p_v / self.dt_speed

    def It_Mean_field_simplified(self, n_rep=None):
        v = self.try_load_property("It_Mean_field_simplified")
        if v is not None:
            return v
        x, y = self.Free_Diff_bis(n_rep=n_rep)[:2]
        # print(self.nori, self.length)
        return x, y * self.nori / self.length * self.p_on * self.p_v / self.dt_speed

    def get_rep_profile(self):
        rep = []
        for il, l in enumerate(self.lengths):
            rep.append(np.zeros(l))
            Nsim = len(self.aRps)
            for sim in range(Nsim):
                rep[il] += np.array(self.aRps[sim][il]) / Nsim
        return rep

    def get_mean_copie(self, time):
        copie = []
        std_copie = []
        rep_t = self.get_times_replication()
        for il, l in enumerate(self.lengths):
            # print(l)
            Nsim = len(self.aRps) - rep_t.tolist().count(-1)
            copie.append(np.ones((Nsim, l)))
            for sim, time_rep in enumerate(rep_t):
                if time_rep != -1:
                    # print("th")
                    copie[il][sim, np.array(self.aRps[sim][il] * self.dte) < time] = 2
                    sim += 1

            std_copie.append(np.std(copie[il], axis=0))
            copie[il] = np.mean(copie[il], axis=0)
        return copie, std_copie

    def Its(self, n_rep=None, recompute=False, cut=0):

        v = self.try_load_property("Its")
        if v is not None:
            return v

        if cut != 0 and recompute is False:
            print("Warning Its does not consider cut")
        elif cut != 0 and recompute is True:
            print("Cut Its considered")

        if recompute:
            NF = self.get_quant("anIts", n_rep=n_rep, cut=cut)[3]
            self.tUNrs = np.sum(np.array(self.aUnrs), axis=1)

            x, _, _, Unr = self.get_quant("tUNrs", n_rep=n_rep)
            Unr[Unr == 0] = np.nan
            y = np.nanmean(NF / Unr, axis=0)
            # Unr[Unr == 0] = 1

            return x, y / self.dt_speed, np.mean(NF, axis=0) / self.dt_speed, np.nanmean(NF, axis=0) / np.nanmean(Unr, axis=0) / self.dt_speed
        else:
            x, y, std, alls = self.get_quant("aIts", n_rep=n_rep)
            # As this are cumulative properties, this scale for one minute
            return x, y / self.dt_speed, std / self.dt_speed, alls / self.dt_speed

    def Ifs(self, n_rep=None, recompute=False, cut=0):
        if recompute == True:
            print("Sorry not the good one implemented")
            return
        if cut != 0 and recompute == False:
            print("Warning Ifs does not consider cut")
        elif cut != 0 and recompute == True:
            print("Cut Ifs considered")

        if recompute:
            self.get_quant("anIts", n_rep=n_rep)
            Nori = self.all + 0
            self.tUNrs = np.sum(np.array(self.aUnrs), axis=1)

            x = self.get_quant("tUNrs", n_rep=n_rep)[0]
            Unr = self.all + 0
            meanurn = np.mean(Unr, axis=0)
            Unr[Unr == 0] = np.nan
            y = np.nanmean(Nori / Unr, axis=0)
            Unr[Unr == np.nan] = 0
            # Unr[Unr == 0] = 1
            return x, y, np.mean(Nori, axis=0), meanurn, Unr
        else:
            return self.get_quant("aIfs", n_rep=n_rep)

    def nIts(self, n_rep=None):
        return self.get_quant("anIts", n_rep=n_rep)

    def MeanIts(self, n_rep=None, cut=0):
        v = self.try_load_property("MeanIts")
        if v is not None:
            return v
        self.tUNrs = np.sum(np.array(self.aUnrs), axis=1)
        x, Nf, std, alls = self.get_quant("anIts", n_rep=n_rep, cut=cut)
        x, Unr, std, allsu = self.get_quant("tUNrs", n_rep=n_rep)
        # allsu[allsu == 0] = np.nan
        print(np.nansum(alls[np.isnan(allsu)]))
        # alls[np.isnan(allsu)] = np.nan
        allsu[np.isnan(allsu)] = 0
        alls[np.isnan(alls)] = 0

        return x, Nf / Unr / self.dt_speed, np.nanmean(alls / allsu, axis=0) / self.dt_speed, np.nanmean(alls, axis=0) / np.nanmean(allsu, axis=0) / self.dt_speed

    def passi_acti(self):
        v = self.try_load_property("passi_acti")
        if v is not None:
            return v

        x, Nori_libre = self.Free_origins()[:2]

        ori_loss = Nori_libre[:-1] - Nori_libre[1:]
        # plot(x[:-1],ori_loss)

        x, activated = self.nIts()[:2]
        # plot(x,activated,label="Activated")

        passivated = ori_loss - activated[:-1]
        # plot(x[:-1],passivated,label="passivated")
        # legend()

        # figure()
        return x[:-1], passivated / activated[:-1]

    def ItsDifferentWay(self, cut=0):
        pass

    def Fds(self, n_rep=None):
        return self.get_quant("aFds", n_rep=n_rep)

    def Free_Diff(self, n_rep=None):
        return self.get_quant("aFree_Diff", n_rep=n_rep)

    def Rps(self, n_rep=None):
        return self.get_quant("aRps", n_rep=n_rep)

    def DNAs(self, n_rep=None):
        v = self.try_load_property("DNAs")
        if v is not None:
            return v
        return self.get_quant("aDNAs", n_rep=n_rep)

    def Free_origins(self, n_rep=None):
        return self.get_quant("aFree_origins", n_rep=n_rep)

    def Free_Diff_bis(self, n_rep=None):
        return self.get_quant("aFree_Diff_bis", n_rep=n_rep)

    def n_activated_oris(self):
        return list(map(len, np.concatenate(self.aFts)))

    def error_DNA_time(self, plot=False, shift=0):

        # https://academic.oup.com/nar/article/42/1/e3/2437422/The-dynamics-of-genome-replication-using-deep
        point = [(4.3714285714285808, 1.0420168067226889), (9.2571428571428562, 1.0126050420168067), (14.40000000000002, 1.0714285714285714), (17.228571428571435, 1.0420168067226889), (19.800000000000015, 0.97058823529411764), (24.428571428571431, 0.96218487394957974), (30.085714285714289, 0.97478991596638642), (32.657142857142873, 1.0714285714285714), (34.71428571428573, 1.1596638655462184), (37.028571428571425, 1.2983193277310923),
                 (39.85714285714284, 1.3277310924369747), (42.428571428571445, 1.3067226890756303), (44.48571428571428, 1.5462184873949578), (46.800000000000026, 1.588235294117647), (49.371428571428581, 1.6470588235294117), (54.771428571428551, 1.672268907563025), (59.914285714285718, 1.8613445378151261), (69.942857142857122, 1.9957983193277311), (79.971428571428589, 1.9495798319327733), (89.742857142857147, 1.8781512605042017)]
        # x_exp,y_exp = zip(*point)

        x, y, std, alls = self.DNAs()
        error = 0
        Np = 0
        for xe, ye in point:
            if xe >= shift:
                i = np.argmin((x - xe + shift)**2)
                # print(x[i],xe)
                error += (ye - y[i])**2
                Np += 1
        if plot:
            return zip(*point)

        return error, Np

    def error_FD_time(self, plot=False, shift=0):

        point = [(250, 0.025), (500, 0.1), (1000, 0.2), (1250, 0.13),
                 (1500, 0.09), (2000, 0.01)]  # Goldar 2008 (/kb)
        point = [(time / 60, value) for time, value in point]

        x, y, std, alls = self.Fds()
        error = 0
        Np = 0
        for xe, ye in point:
            if xe >= shift:
                i = np.argmin((x - xe + shift)**2)
                # print(x[i],xe)
                error += (ye - y[i])**2
                Np += 1
        if plot:
            return zip(*point)

        return error, Np

    def error_firing_time(self, plot=False, specie="yeast", coarse=1, arach=False, smooth=1):

        # Universal Temporal Prrofile of Replication Origin (Goldar)
        if not specie in ["yeast", "xenope"]:
            raise
        point = [(5, 0.01), (13, 0.02), (16, 0.04), (20, 0.07), (25, 0.02),
                 (30, 0.01)] + [(i, 0) for i in range(31, 70, 2)]  # xenoput
        unity = 1  # we want it by minutes
        point = [(time, value * unity) for time, value in point]
        if specie == "yeast":
            point = [11.104005791505799, 0.00018581081081081065,
                     12.066008316008308, 0.00020270270270270323,
                     13.165837540837543, 0.00023648648648648667,
                     13.990477427977439, 0.0002533783783783784,
                     15.0921629046629, 0.0003547297297297296,
                     16.05787793287793, 0.0005067567567567568,
                     17.161883724383713, 0.0006925675675675674,
                     18.127134689634687, 0.0008277027027027029,
                     19.092849717849717, 0.0009797297297297301,
                     20.19592738342739, 0.0011317567567567573,
                     21.159786159786165, 0.001216216216216216,
                     22.1227168102168, 0.001266891891891892,
                     23.22393822393822, 0.0013513513513513514,
                     24.191509504009503, 0.001570945945945946,
                     25.298763736263723, 0.001875,
                     26.407410157410155, 0.0022297297297297295,
                     27.233442233442233, 0.0022972972972972973,
                     28.46970596970597, 0.0022972972972972973,
                     29.431244431244423, 0.0022972972972972973,
                     30.402528215028198, 0.0026520270270270273,
                     31.514887139887136, 0.0031418918918918915,
                     32.35437704187704, 0.003699324324324324,
                     33.59156890406891, 0.003733108108108108,
                     34.55125111375111, 0.0036655405405405404,
                     35.50907707157708, 0.003530405405405405,
                     36.614475051975035, 0.0037668918918918916,
                     37.723121473121466, 0.004121621621621621,
                     38.69208494208493, 0.004391891891891891,
                     39.65640778140778, 0.004493243243243243,
                     40.747419809919805, 0.004206081081081081,
                     41.696892634392626, 0.0037668918918918916,
                     42.666320166320176, 0.004054054054054054,
                     43.775894713394706, 0.004442567567567567,
                     44.73279254529254, 0.004273648648648648,
                     45.82380457380458, 0.003986486486486486,
                     46.62338506088507, 0.003091216216216216,
                     47.83180501930502, 0.0020777027027027027,
                     48.78591847341846, 0.0018074324324324326,
                     49.72425378675379, 0.0009628378378378375,
                     50.65934065934067, 0,
                     51.75824175824175, 0,
                     52.85760692010692, 0.000016891891891892587,
                     53.81914538164537, 0.000016891891891892587,
                     54.780219780219795, 0,
                     56.15384615384616, 0,
                     57.11538461538461, 0,
                     57.93956043956044, 0]
            point = np.array(point)
            point = point.reshape(-1, 2)

        if arach:
            if specie == "yeast":
                print(point.shape)
                point = pd.read_csv("../../data/external/I2T_yeast.txt", sep="  ", header=None)
                point[0] += 10
                point = np.array(point)
                print(point.shape)

            if specie == "xenope":
                A = pd.read_csv("../../data/external/I2T_exenope.txt", sep="  ", header=None)

                A[0] = (A[0] - 20 * 60) / 60
                A[1] = pd.rolling_mean(A[1], window=smooth) * 10
                point = np.array(A)

        x, y = self.Its()[:2]
        error = 0
        Np = 0
        shift = 0
        for xe, ye in point:
            if xe >= shift:
                i = np.argmin((x - xe + shift)**2)
                # print(x[i],xe)
                error += (ye - y[i] / coarse)**2
                Np += 1
        if plot:
            return zip(*point)
        return error, Np

    def xenope_prof(self, profile=True, which="mean"):
        import matplotlib.pyplot as plt

        chro = 0
        coarse = 1000
        if profile:
            if which == "mean":
                Prof = self.get_rep_profile()[chro]
                x = np.arange(len(Prof)) * coarse / 1000.
                y = Prof * self.dte
                plt.plot(x, Prof * self.dte, label="Simulated")
                plt.xlim(-10, x[-1] + 10)
            else:
                for sim in which:
                    x = np.arange(len(self.aRps[sim][chro])) * coarse / 1000.
                    plt.plot(x, self.aRps[sim][chro] * self.dte)
                top = self.aRps[sim][chro]
                plt.xlim(-10, x[-1] + 10)
        else:
            k = list(times.keys())
            k.sort()
            for ikk, kk in enumerate(k):
                if ikk == 0:
                    mean_C = mean_copie[kk][chro]
                else:
                    mean_C += mean_copie[kk][chro]
            x = np.arange(len(mean_C)) * coarse / 1000.
            plt.plot(np.arange(len(mean_C)) * coarse / 1000., mean_C / len(k))
            plt.xlim(-10, x[-1] + 10)
        return x, y

    def whole_genome_timing(self, coarse=5000, figsize=(12, 12), plot=True,
                            default_rep="../../data/external/time-coordinate.pick",
                            experiment=True, profile=False, which="mean", fig=None,
                            warning=True, ori=True, shift=0, N_chrom=range(16), strength_ori=None):

        import matplotlib.pyplot as plt

        with open(default_rep, "rb") as f:
            times, coordinate = cPickle.load(f)

        times.keys()
        time_p = list(times.keys())
        time_p.sort()
        dna = []
        for t in time_p:
            dna.append(np.concatenate(times[t]).mean())
        # plot(time_p, dna)

        result = {"chr": [], "start": [], "end": [], "mean_copie_exp": [], "mean_copie_simu": []}
        # f = figure(figsize=(20,20))
        if fig is None:
            f = plt.figure(figsize=figsize)
        elif fig == "other":
            pass
        else:
            f = fig
        mean_copie = {}
        if not profile:
            k = list(times.keys())
            k.sort()
            for ikk, kk in enumerate(k):

                mean_copie[kk] = self.get_mean_copie(max(0, int(kk) - shift))[0]
                # print(mean_copie[kk],len(mean_copie[kk][0]) )
                # print(len( mean_copie[kk]))

        if profile:
            max_t = self.get_times_replication()
            if which == "mean":

                max_t = max(max_t)
            else:

                max_t = max(np.array(max_t)[which])
            if max_t == -1:
                max_t = np.max(self.get_times_replication(finished=False))
        extra = [0, 0, 0, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6]
        position = [0, 1, 2, 0, 0, 1, 2, 0, 1, 2, 0, 1, 0, 1, 0, 1]
        s = 0.03
        sh = 0.04
        height = 1 / (7 + 1) - sh

        margin_right = 0.02

        for chro in N_chrom:
            # ax = f.add_subplot(4,4,chro + 1)
            # ax = f.add_subplot(gs[chro])

            # print([xstart,ystart,w,h])
            if fig != "other":
                column = extra[chro]
                tot = extra.count(column)
                p = position[chro]

                row_lengths = [l for l, i in zip(self.lengths, extra) if column == i]
                crow_length = [0] + np.cumsum(row_lengths).tolist()

                xstart = (p + 1) * s + (1 - margin_right - tot * s) * \
                    crow_length[p] / (sum(row_lengths))
                ystart = 1 - (column + 1) * (height + sh)
                w = (1 - margin_right - tot * s) * row_lengths[p] / (sum(row_lengths))
                h = height
                f.add_axes([xstart, ystart, w, h])

            # chro = 3
            if profile:
                if which == "mean":
                    Prof = self.get_rep_profile()[chro]
                    x = np.arange(len(Prof)) * coarse / 1000.
                    plt.plot(x, Prof * self.dte, label="Simulated")
                    plt.xlim(-10, x[-1] + 10)
                else:
                    for sim in which:
                        x = np.arange(len(self.aRps[sim][chro])) * coarse / 1000.
                        plt.plot(x, self.aRps[sim][chro] * self.dte)
                    top = self.aRps[sim][chro]
                    plt.xlim(-10, x[-1] + 10)
            else:
                k = list(times.keys())
                k.sort()
                for ikk, kk in enumerate(k):
                    if ikk == 0:
                        mean_C = mean_copie[kk][chro]
                    else:
                        mean_C += mean_copie[kk][chro]
                x = np.arange(len(mean_C)) * coarse / 1000.
                plt.plot(np.arange(len(mean_C)) * coarse / 1000., mean_C / len(k))
                plt.xlim(-10, x[-1] + 10)

                top = mean_C / len(k)
            if ori:
                if strength_ori is not None:
                    st = strength_ori[chro]
                else:
                    st = [1] * len(self.l_ori[chro])
                for x, s in zip(self.l_ori[chro], st):
                    # print(np.array(top)[~np.equal(top, None)])

                    mini = min(np.array(top)[~np.equal(top, None)])
                    maxi = max(np.array(top)[~np.equal(top, None)])
                    mini = 1
                    maxi = 2
                    # print(mini, maxi)
                    col = "k"
                    if s != 1:
                        col = {"Confirmed": "r", "Likely": "g", "Dubious": "b"}[s]
                    plt.plot([x * coarse / 1000., x * coarse / 1000],
                             [mini, maxi], "--", color=col, linewidth=1)

            def get_rep_prof(times, coordinate, ch, profile=True):
                k = list(times.keys())
                k.sort()

                # To get all the coordinates
                m = []
                for kk in k:
                    m = list(set(coordinate[kk][ch] + m))
                m.sort()
                # print(len(m))

                rep = np.zeros(len(m))  # + 70
                norm = np.zeros(len(m))
                for ilocus, locus in enumerate(m):
                    # print(locus)
                    for kk in k[:: -1]:
                        if locus in coordinate[kk][ch]:
                            i = list(coordinate[kk][ch]).index(locus)
                            if profile:
                                if times[kk][ch][i] > 1.5:
                                    rep[ilocus] = min(int(kk), 70)

                            else:
                                # Mean replication value
                                rep[ilocus] += times[kk][ch][i]
                                norm[ilocus] += 1
                norm[norm == 0] = 1
                if profile:
                    rep[rep == 0] = 70
                # print(times[kk][ch])
                return m, rep / norm
            if experiment:
                locci, p = get_rep_prof(times, coordinate, chro, profile=profile)

                # m = lengths[chro] / len(p)
                # plot(np.arange(len(p)) * m,p)
                if not profile:
                    for loc, copie in zip(locci, p):
                        result["chr"].append(chro + 1)

                        result["start"].append(loc)
                        result["end"].append(loc)
                        result["mean_copie_exp"].append(copie)
                        try:
                            result["mean_copie_simu"].append(top[int(loc / coarse)])
                        except IndexError:
                            if warning:
                                print("out of bounds")
                            result["mean_copie_simu"].append(top[-1])

                plt.plot(np.array(locci) / 1000., p, "-", label="simulated")
            if profile:
                plt.ylim(max_t, 0)

            else:
                plt.ylim(1., 2.)
            if extra[chro] == 6:
                plt.xlabel("Genomic position (kb)")
            if position[chro] == 0:
                if profile:
                    plt.ylabel("rep time (min)")
                else:
                    plt.ylabel("gene copy number")


class ensembleSimAnalysis(ensembleSim):

    def __init__(self, json_file, hdf5_file):
        with open(json_file, "r") as f:
            self.parameters = json.load(f)

        sub_sample_ori = self.parameters.get("sub_sample_ori", None)
        if sub_sample_ori:
            self.parameters.pop("sub_sample_ori")
            l_ori = [list(range(int(self.parameters["lengths"][0] * sub_sample_ori)))]
            lengths = self.parameters["lengths"]
        else:
            extra = "../../"
            if type(self.parameters["lengths"]) == str:
                lengths, _ = load_lengths_and_centro(
                    extra + self.parameters["lengths"], self.parameters["coarse"])
                self.parameters["lengths"] = lengths

            if type(self.parameters["Nori"]) == str and self.parameters["Nori"] != "xenope":

                d = {"C": "Confirmed", "L": "Likely", "D": "Dubious"}
                ot = []
                for o in self.parameters["ori_type"]:
                    ot.append(d[o])

                l_ori = load_ori_position(extra + self.parameters["Nori"],
                                          ot,
                                          self.parameters["lengths"],
                                          self.parameters["coarse"])

        ensembleSim.__init__(self, Nsim=self.parameters["Nsim"],
                             Nori=None, Ndiff=self.parameters["Ndiff"],
                             lengths=self.parameters["lengths"],
                             p_on=self.parameters["p_on"],
                             p_v=self.parameters["p_v"],
                             dt_speed=self.parameters["dt_speed"],
                             fork_speed=self.parameters["fork_speed"],
                             p_off=None, only_one=True, l_ori=l_ori)
        self.hdf5_file = hdf5_file
