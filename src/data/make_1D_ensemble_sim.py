
import sys
sys.path.append("./")
from replication.ensembleSim import ensembleSim
from replication.simulate import load_parameters
from replication.tools import load_ori_position, load_lengths_and_centro
import os
import json
import _pickle as cPickle
import numpy as np
import time
import os
import errno


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

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
    make_sure_path_exists(parameters["data_folder"])
    with open(os.path.join(parameters["data_folder"], "params.json"), "w") as f:
        s = json.dumps(parameters)
        f.write(s)

    # ensembleSim(Nsim, Nori, Ndiff, lengths, p_on, p_off, only_one,
    # all_same_ori=False, l_ori=[], cut=10)
    if type(parameters["lengths"]) == str:
        lengths, _ = load_lengths_and_centro(parameters["lengths"], parameters["coarse"])
        parameters["lengths"] = lengths

    if type(parameters["Nori"]) == str and parameters["Nori"] != "xenope":

        d = {"C": "Confirmed", "L": "Likely", "D": "Dubious"}
        ot = []
        for o in parameters["ori_type"]:
            ot.append(d[o])

        coarsed = parameters.get("coarsed", True)
        l_ori = load_ori_position(parameters["Nori"],
                                  ot,
                                  parameters["lengths"],
                                  parameters["coarse"], coarsed=coarsed)

        if not coarsed:
            parameters["positions"] = [] + l_ori
            l_ori = [list(range(len(ip_ori))) for ip_ori in parameters["positions"]]
        if "coarsed" in parameters:
            c = parameters.pop("coarsed")
            print(c, sum(list(map(len, l_ori))))

    if "coarse" in parameters:
        parameters.pop("coarse")

    if parameters["Nori"] == "xenope":
        sub_sample_ori = parameters.pop("sub_sample_ori")
        l_ori = [list(range(int(parameters["lengths"][0] * sub_sample_ori)))]

        positions = [[]]
        if parameters["homogeneous"]:
            interval = parameters["lengths"][0] / 1.0 / len((l_ori[0]))
            for i in range(len((l_ori[0]))):
                positions[-1].append(i * interval + interval * np.random.uniform())
        else:
            positions = []
            for sim in range(parameters["Nsim"]):
                spositions = [[]]
                for i in range(len((l_ori[0]))):
                    spositions[0].append(parameters["lengths"][0] * np.random.uniform())
                spositions[0].sort()
                positions.append(spositions)
        parameters.pop("homogeneous")
        parameters["positions"] = positions

    parameters.pop("filename")
    data_folder = parameters.pop("data_folder")
    parameters.pop("ori_type")

    parameters["Nori"] = l_ori

    correlation = parameters.get("correlation", True)
    if "correlation" in parameters:
        parameters.pop("correlation")

    if "ramp_time" in parameters:
        if parameters.get("ramp_type", "linear") == "linear":

            parameters["max_ramp"] = parameters["Ndiff"]
            # if parameters["Nori"] == "xenope":
            parameters["ramp"] = parameters["max_ramp"] / parameters.pop("ramp_time")

        if parameters.get("ramp_type") == "exp":
            print(parameters)
            parameters["ramp"] = parameters.pop("ramp_time")
            parameters["max_ramp"] = parameters["Ndiff"]

    if "ramp_type" in parameters:
        if parameters["ramp_type"] == "pulse":
            parameters["ramp"] = False
    # else:

    E = ensembleSim(**parameters)
    E.run_all(20000, correlation=correlation)
    with open(os.path.join(data_folder, "ensembleSim.pick"), "wb") as f:
        cPickle.dump(E, f)
