
import sys
sys.path.append("./")
import os
import errno
import glob
import _pickle as cPickle
import sys
sys.path.append("../../src/data")

from replication.tools import load_3D_simus


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

if __name__ == "__main__":
    simu = sys.argv[1]

    threeD = False
    if len(sys.argv) == 3:
        threeD = bool(sys.argv[2])
    # print(sys.argv)

    root = "data/raw"
    two = True

    if not threeD:
        with open(root + "/" + simu + "/" + "ensembleSim.pick", "rb") as f:
            tmp_simu = cPickle.load(f)
    else:
        ni = len(glob.glob(root + "/" + simu + "/traj*"))
        tmp_simu, lengths, parameters1 = load_3D_simus(root + "/" + simu + "/traj", n=ni, orip=True)
        tmp_simu.Nsim = len(tmp_simu.aIts)

        tmp_simu.add_precomputed("Mean_replication_time", root + "/" + simu + "/" + "analysis.hdf5")
        tmp_simu.add_precomputed("get_rep_profile", root + "/" + simu + "/" + "analysis.hdf5")
        tmp_simu.add_precomputed("n3Dsim", root + "/" + simu + "/" + "analysis.hdf5")

    tmp_simu.add_precomputed("passi_acti", root + "/" + simu + "/" + "analysis.hdf5")
    tmp_simu.add_precomputed("passi", root + "/" + simu + "/" + "analysis.hdf5")
    tmp_simu.add_precomputed("acti", root + "/" + simu + "/" + "analysis.hdf5")

    tmp_simu.add_precomputed("Its", root + "/" + simu + "/" + "analysis.hdf5", two=two)
    tmp_simu.add_precomputed("nIts", root + "/" + simu + "/" + "analysis.hdf5", two=two)
    tmp_simu.add_precomputed("rho_ori", root + "/" + simu + "/" + "analysis.hdf5", two=two)
    tmp_simu.add_precomputed("Free_Diff", root + "/" + simu + "/" + "analysis.hdf5", two=two)
    tmp_simu.add_precomputed("Free_Diff_bis", root + "/" +
                             simu + "/" + "analysis.hdf5", two=two)

    tmp_simu.add_precomputed("Fds", root + "/" + simu + "/" + "analysis.hdf5", two=two)

    tmp_simu.add_precomputed("get_times_replication", root + "/" + simu + "/" + "analysis.hdf5")
    tmp_simu.add_precomputed("get_dist_between_activated_origins",
                             root + "/" + simu + "/" + "analysis.hdf5")
    tmp_simu.add_precomputed("MeanIts", root + "/" + simu + "/" + "analysis.hdf5", two=two)
    tmp_simu.add_precomputed("It_Mean_field_origins", root + "/" + simu + "/" + "analysis.hdf5")
    tmp_simu.add_precomputed("It_Mean_field_simplified", root + "/" +
                             simu + "/" + "analysis.hdf5")

    tmp_simu.add_precomputed("DNAs", root + "/" + simu + "/" + "analysis.hdf5", two=two)
