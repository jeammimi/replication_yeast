import sys
sys.path.append("./")
from replication.ensembleSim import ensembleSim
from replication.simulate import load_parameters
from replication.tools import load_ori_position, load_lengths_and_centro
import os
import json
import _pickle as cPickle
from skopt import gp_minimize
from skopt import dump
import numpy as np


def latin(n, ranges, save=False):
    """
    Build latin hypercube.
    Parameters
    ----------
    n : int
        Number of points.
    d : int
        Size of space.
    Returns
    -------
    pts : ndarray
        Array of points uniformly placed in d-dimensional unit cube.
    """
    # starting with diagonal shape
    d = len(ranges)
    pts = np.ones((n, d))

    for i in range(n):
        pts[i] = pts[i] * i / (n - 1.)

    # spread function
    def spread(p):
        s = 0.
        for i in range(n):
            for j in range(n):
                if i > j:
                    s = s + 1. / np.linalg.norm(np.subtract(p[i], p[j]))
        return s

    # minimizing spread function by shuffling
    currminspread = spread(pts)
    if save:
        Save = [pts]

    for m in range(1000):

        p1 = np.random.randint(n)
        p2 = np.random.randint(n)
        k = np.random.randint(d)

        newpts = np.copy(pts)
        newpts[p1, k], newpts[p2, k] = newpts[p2, k], newpts[p1, k]
        newspread = spread(newpts)

        if newspread < currminspread:
            pts = np.copy(newpts)
            currminspread = newspread
            if save:
                Save.append(np.copy(newpts))

    for ir, r in enumerate(ranges):
        print(r)
        pts[::, ir] = (r[1] - r[0]) * pts[::, ir] + r[0]

    if save:
        return pts, Save

    return pts

if __name__ == "__main__":
    param_file = sys.argv[1]
    parameters = load_parameters(param_file)
    # print(sys.argv)
    if "sumatra_label" in parameters:
        parameters["data_folder"] = os.path.join(parameters["data_folder"],
                                                 parameters["sumatra_label"])
        parameters.pop("sumatra_label")
    else:
        print("no extra label")

    parameters["data_folder"] = os.path.join(parameters["data_folder"], "")
    parameters["filename"] = param_file

    print(parameters["data_folder"])
    with open(os.path.join(parameters["data_folder"], "params.json"), "w") as f:
        s = json.dumps(parameters)
        f.write(s)

    # ensembleSim(Nsim, Nori, Ndiff, lengths, p_on, p_off, only_one,
    # all_same_ori=False, l_ori=[], cut=10)
    if type(parameters["lengths"]) == str:
        lengths, _ = load_lengths_and_centro(parameters["lengths"], parameters["coarse"])
        parameters["lengths"] = lengths

    if type(parameters["Nori"]) == str and parameters["Nori"] != "xenope":
        l_ori = load_ori_position(parameters["Nori"],
                                  parameters["ori_type"],
                                  parameters["lengths"],
                                  parameters["coarse"])

    if parameters["Nori"] == "xenope":
        l_ori = [list(range(parameters["lengths"][0]))]

    parameters.pop("filename")
    data_folder = parameters.pop("data_folder")
    parameters.pop("ori_type")
    c = parameters.pop("coarse")

    parameters["Nori"] = l_ori

    def error(x, returnv=False, c=c):
        if len(x) == 2:
            only_one = True
        else:
            only_one = False
        if only_one:
            Ndiff = x[0]
            p_on = x[1]
            p_off = 0.2
        else:
            Ndiff = x[0]
            p_on = x[1]
            p_off = x[2]
        # print(Ndiff,p_on)
            # print(Ndiff,p_on)
        E = ensembleSim(parameters["Nsim"], l_ori, Ndiff, parameters["lengths"],
                        p_on=p_on, p_off=p_off, only_one=only_one,
                        dt_speed=parameters["dt_speed"],
                        fork_speed=parameters["fork_speed"],
                        gindin=parameters["gindin"],
                        p_v=parameters["p_v"],
                        random=True,
                        all_same_ori=True)
        E.run_all(40000)

        if returnv:
            return E
        else:
            if parameters["Nori"] == "xenope":
                c = 1
            else:
                c = c / 1000
            # print(c)
            return getattr(E, parameters["optimisation"])(coarse=c)[0]

    if parameters["only_one"]:
        start = [parameters["rNdiff"], parameters["rp_on"]]
    else:
        start = [parameters["rNdiff"], parameters["rp_on"], parameters["rp_off"]]
    x0 = latin(10, start)
    print(x0)
    res = gp_minimize(error, start, n_jobs=1, n_calls=100, n_random_starts=0, x0=x0)

    with open(os.path.join(data_folder, "ensembleSim.pick"), "wb") as f:
        cPickle.dump(error(res["x"], returnv=True), f)

    dump(res, os.path.join(data_folder, "optimisation.pkl"))
