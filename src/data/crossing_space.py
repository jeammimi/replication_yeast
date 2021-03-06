from replication.simulate import create_initial_configuration
from replication.simulate import load_parameters, minimize
from hoomd import init, group, md, deprecated, dump, analyze
import hoomd
import numpy as np
import os
import sys
import json

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
    with open(os.path.join(parameters["data_folder"], "params.json"), "w") as f:
        s = json.dumps(parameters)
        f.write(s)

    hoomd.context.initialize()  # "--mode=cpu ")

    # print(type(traj["p_origins"]) == list)
    # if hoomd.comm.get_rank() == 0:
    traj = parameters
    snapshot, _, tag_spb, bond_list, plist, Cp, lP = create_initial_configuration(traj)


    R = traj["R"]
    data_folder = traj["data_folder"]
    dcd_period = traj["dcd_period"]
    control = traj["control"]


    plist = ["A", "B"]
    bond_list = ["A-A"]
    snapshot.particles.types = plist
    snapshot.bonds.types = bond_list
    for p in range(len(snapshot.particles.typeid)):
        snapshot.particles.typeid[p] = np.random.randint(2)
    for p in range(len(snapshot.bonds.typeid)):
        snapshot.bonds.typeid[p] = 0

    system = init.read_snapshot(snapshot)

    xml = deprecated.dump.xml(
        filename=data_folder +
        "atoms.hoomdxml",
        period=None,
        group=group.all(),
        vis=True)

    logger = analyze.log(
        filename=data_folder +
        'mylog.log',
        period=1000,
        quantities=[
            'temperature',
            'potential_energy',
            'bond_harmonic_energy',
            'external_wall_lj_energy',
            "pair_table_energy",
            'kinetic_energy',
            'volume',
            'pressure'],
        overwrite=True)
    # xml.disable()

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
    if not control:
        table.pair_coeff.set("A", "B",
                             func=cos_soft, rmin=0, rmax=r_cut,
                             coeff=dict(epsilon=0, sigma=1.0))

    sphere = md.wall.group()
    sphere.add_sphere(r=R, origin=(0.0, 0.0, 0.0), inside=True)
    # lj much more slower (at least in thu minimisation)
    wall_force_slj = md.wall.lj(sphere, r_cut=1.12)
    wall_force_slj.force_coeff.set(plist, epsilon=1.0, sigma=1.0,
                                   r_cut=1.12, mode="shift")
    # hoomd.comm.barrier()
    microtubule_length = None
    Spb_g = None
    all_move = group.all()
    minimize(traj, all_move, system, snapshot, Spb_g, Cp, microtubule_length)

    sim_dt = traj["sim_dt"]
    seed = traj["seed"]

    md.integrate.mode_standard(dt=sim_dt)
    method = md.integrate.langevin(group=all_move, kT=1, seed=seed)
    snp = system  # .take_snapshot()


    md.integrate.mode_standard(dt=sim_dt / 4)
    hoomd.run(100)
    md.integrate.mode_standard(dt=sim_dt / 2)
    hoomd.run(100)

    dcd = dump.dcd(filename=data_folder + 'poly.dcd',
                   period=dcd_period, overwrite=True)
    hoomd.run(1000000)

    dcd.disable()
    logger.disable()
