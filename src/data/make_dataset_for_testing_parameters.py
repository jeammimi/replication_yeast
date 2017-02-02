import sys
sys.path.append("./")
from replication.simulate import simulate, load_parameters
import os
import json

if __name__ == "__main__":
    param_file = sys.argv[1]
    parameters = load_parameters(param_file)
    print(sys.argv)
    if len(sys.argv) >= 3:
        parameters["visu"] = True
        if "sumatra_label" in parameters:
            parameters.pop("sumatra_label")

    if "sumatra_label" in parameters:
        parameters["data_folder"] = os.path.join(parameters["data_folder"],
                                                 parameters["sumatra_label"])
        parameters["data_folder"] = os.path.join(parameters["data_folder"], "")
    else:
        print("no extra label")

    parameters["filename"] = param_file
    with open("params.json", "w") as f:
        s = json.dumps(parameters)
        f.write(s)

    simulate(parameters)
