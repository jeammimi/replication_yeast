
import sys
sys.path.append("./")
import os
import errno

import _pickle as cPickle
import sys
sys.path.append("../../src/data")

from replication.ensembleSim import ensembleSim


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

if __name__ == "__main__":
    simu = sys.argv[1]
    # print(sys.argv)

    root = "data/raw"
    two = True
    with open(root + "/" + simu + "/" + "ensembleSim.pick", "rb") as f:

        tmp_simu = cPickle.load(f)
        tmp_simu.add_precomputed("passi_acti", root + "/" + simu + "/" + "analysis.hdf5")
        tmp_simu.add_precomputed("Its", root + "/" + simu + "/" + "analysis.hdf5", two=two)
        tmp_simu.add_precomputed("get_times_replication", root + "/" + simu + "/" + "analysis.hdf5")
        tmp_simu.add_precomputed("get_dist_between_activated_origins",
                                 root + "/" + simu + "/" + "analysis.hdf5")
        tmp_simu.add_precomputed("MeanIts", root + "/" + simu + "/" + "analysis.hdf5", two=two)
        tmp_simu.add_precomputed("It_Mean_field_origins", root + "/" + simu + "/" + "analysis.hdf5")
        tmp_simu.add_precomputed("It_Mean_field_simplified", root + "/" +
                                 simu + "/" + "analysis.hdf5")

        tmp_simu.add_precomputed("DNAs", root + "/" + simu + "/" + "analysis.hdf5", two=two)
