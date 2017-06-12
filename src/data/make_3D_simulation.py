import sys
sys.path.append("./")
from replication.simulate import simulate, load_parameters
import os
import json
import copy


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
    if len(sys.argv) >= 3:
        parameters["visu"] = True
        if "sumatra_label" in parameters:
            parameters.pop("sumatra_label")

    if "sumatra_label" in parameters:
        parameters["data_folder"] = os.path.join(parameters["data_folder"],
                                                 parameters["sumatra_label"])
    else:
        print("no extra label")

    parameters["data_folder"] = os.path.join(parameters["data_folder"], "")
    parameters["filename"] = param_file

    print(parameters["data_folder"])
    make_sure_path_exists(parameters["data_folder"])
    with open(os.path.join(parameters["data_folder"], "params.json"), "w") as f:
        s = json.dumps(parameters)
        f.write(s)

    original = copy.deepcopy(parameters["visu"])
    parameters["visu"] = True
    simulate(parameters)  # generate files for visualisation
    parameters["visu"] = original
    simulate(parameters)
