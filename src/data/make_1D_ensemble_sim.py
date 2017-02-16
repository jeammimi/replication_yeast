
import sys
sys.path.append("./")
from replication.ensembleSim import ensembleSim
from replication.simulate import load_parameters
import os
import json
import _pickle as cPickle
import pandas


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
    if type(parameters["Nori"]) == str:
        Oris = pandas.read_csv(parameters["Nori"], comment="#")

        l_ori = [[] for i in range(16)]

        # Filters ori:
        for ch, p, status in zip(Oris["chr"], Oris["start"], Oris["status"]):
            if status in parameters["ori_type"]:  # ,"Likely"]:
                l_ori[ch - 1].append(int(p / 1000))

        # Then remove duplicate (/kb ) and outside of boundaries

        tot = 0
        for i in range(len(parameters["lengths"])):
            isize = len(l_ori[i])
            l_ori[i] = list(set(l_ori[i]))
            l_ori[i].sort()
            print(isize, len(l_ori[i]))
            while not (max(l_ori[i]) < parameters["lengths"][i]):
                l_ori[i].remove(max(l_ori[i]))
            # print(max(l_ori[i]),len_chrom[i]*5

            tot += len(l_ori[i])
            assert(max(l_ori[i]) < parameters["lengths"][i])

    parameters.pop("filename")
    data_folder = parameters.pop("data_folder")

    parameters["Nori"] = l_ori
    E = ensembleSim(**parameters)
    E.run_all(200)
    with open(os.path.join(data_folder, "ensembleSim.pick"), "wb") as f:
        cPickle.dump(E, f)
