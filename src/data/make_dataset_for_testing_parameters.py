import sys
sys.path.append("./")
from replication.simulate import simulate, load_parameters
import os

if __name__ == "__main__":
    param_file = sys.argv[1]
    parameters = load_parameters(param_file)
    parameters["data_folder"] = os.path.join(parameters["data_folder"],
                                             parameters["sumatra_label"])
