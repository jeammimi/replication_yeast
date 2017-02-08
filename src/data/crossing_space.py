from replication.simulate import create_initial_configuration, force_field
from replication.simulate import load_parameters, minimize
from hoomd import init, group, md
import hoomd
import numpy as np

if __name__ == "__main__":

    traj = load_parameters("replication/cross.json")
    traj["data_folder"] = "/tmp/"
    traj["p_origins"] = [[10] for i in range(16)]

    hoomd.context.initialize()  # "--mode=cpu ")

    # print(type(traj["p_origins"]) == list)
    # if hoomd.comm.get_rank() == 0:

    snapshot, _, tag_spb, bond_list, plist, Cp, lP = create_initial_configuration(traj)

    plist = ["A", "B"]
    bond_list = ["A-A"]
    for p in range(len(snapshot.particles.typeid)):
        snapshot.particles.typeid[p] = np.random.randint(2)
    for p in range(len(snapshot.bonds.typeid)):
        snapshot.bonds.typeid[p] = 0

    system = init.read_snapshot(snapshot)

    # Force_field:

    harmonic = md.bond.harmonic()
    harmonic.bond_coeff.set(bond_list, k=20.0, r0=1)

    def cos_soft(r, rmin, rmax, epsilon, sigma):

        V = epsilon * (1 + np.cos(r * 3.1415 / (rmax)))
        F = epsilon * 3.1415 / (rmax) * np.sin(r * 3.1415 / (rmax))

        return (V, F)


    nl = md.nlist.tree(r_buff=0.4, check_period=1)
    # nl = md.nlist.stencil(r_buff=0.4, check_period=1)
    # nl = md.nlist.cell(r_buff=0.4, check_period=1)

    r_cut = 1.5
    epsilon = 6.5

    table = md.pair.table(width=1000, nlist=nl)
    table.pair_coeff.set(plist, plist,
                         func=cos_soft, rmin=0, rmax=r_cut,
                         coeff=dict(epsilon=epsilon, sigma=1.0))
    table.pair_coeff.set("A", "B",
                         func=cos_soft, rmin=0, rmax=r_cut,
                         coeff=dict(epsilon=0, sigma=1.0))

    # hoomd.comm.barrier()
    microtubule_length = None
    Spb_g = None
    all_move = group.all()
    minimize(traj, all_move, system, snapshot, Spb_g, Cp, microtubule_length)
