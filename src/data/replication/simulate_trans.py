# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 16:36:10 2017

@author: jarbona
"""

import hoomd
from hoomd import data, init, md, group, dump, deprecated, analyze
import numpy as np
import scipy.linalg as linalg
from scipy.spatial.distance import cdist
from .PMotion import Polymer, Diffusing
import _pickle as cPickle
from .createPoly import create_init_conf_yeast
from replication.tools import load_ori_position, load_lengths_and_centro
from replication.tools import load_parameters
import time
import json


def create_initial_configuration(traj):

    len_chrom = traj["len_chrom"]
    Cent = traj["Cent"]
    p_ribo = traj["p_ribo"]
    R = traj["R"]
    micron = traj["micron"]

    # Diffusing elements
    N_diffu = traj["N_diffu"]
    p_origins = traj["p_origins"]

    if type(len_chrom) != list:
        len_chrom, _ = load_lengths_and_centro(len_chrom, traj["coarse"])

    if type(Cent) != list:
        _, Cent = load_lengths_and_centro(Cent, traj["coarse"])

    if type(p_origins) != list:
        p_origins = load_ori_position(traj["p_origins"],
                                      traj["ori_type"],
                                      len_chrom,
                                      traj["coarse"])
    p_ribo = [[int(position) // int(traj["coarse"]), length] for position, length in p_ribo]

    two_types = traj.get("two_types", False)
    p_second = traj.get("p_second", [])
    dstrength = traj.get("dstrength", 0)
    strengths = None
    if p_second != []:

        # Assign delta_strength
        if dstrength != 0:
            strengths = []
            for bands, pos in zip(p_second, p_origins):
                strengths.append([])
                for p in pos:
                    found = False
                    for Intervals in bands:
                        if Intervals[0] < p < Intervals[1]:
                            strengths[-1].append(dstrength)
                            found = True
                            break
                    if not found:
                        strengths[-1].append(1)

        ps = []
        for ch in p_second:
            ps.append([])
            for p1, p2 in ch:
                ps[-1] += range(p1, p2)
        p_second = ps

    print("strengths", strengths)
    # Yeast case
    spb = traj["spb"]
    nucleole = traj["nucleole"]
    telomere = traj["telomere"]
    microtubule_length = traj["microtubule_length"] * micron
    special_start = traj["special_start"]
    visu = traj["visu"]
    dump_hic = traj["dump_hic"]

    # Scenari
    diff_bind_when_free = traj["diff_bind_when_free"]

    # Simulation parameters

    Np = len(len_chrom)
    assert(len(len_chrom) == len(Cent) == len(p_ribo))
    if special_start:
        Sim = create_init_conf_yeast(
            len_chrom=len_chrom,
            dist_centro=Cent,
            p_ribo=p_ribo,
            Radius=R - 1,
            Mt=microtubule_length)
    else:
        Sim = []

    spbp = 0 if not spb else 1

    Total_particle = sum(len_chrom) + N_diffu * 2 + spbp
    list_nuc = [list(range(start, start + size)) if size != 0 else []
                for start, size in p_ribo]
    # print(list_nuc)
    # exit()

    snapshot = data.make_snapshot(
        N=Total_particle, box=data.boxdim(
            L=2 * R), bond_types=['polymer'])

    spbb = Np if spb else 0

    if visu:
        spbb = 0

    bond_diffu = 0
    if diff_bind_when_free:
        bond_diffu = N_diffu

    snapshot.bonds.resize(sum(len_chrom) - len(len_chrom) + bond_diffu + spbb)

    bond_list = ['Mono_Mono', 'Diff_Diff', 'Mono_Diff']

    if spb:
        bond_list += ["Spb_Cen"]
    if nucleole:
        bond_list += ["Mono_Nuc", "Nuc_Nuc"]

    snapshot.bonds.types = bond_list

    plist = ['Mono', 'Ori', 'Diff', 'S_Diff', 'F_Diff', "I_Diff"]
    if two_types:
        plist.append("Mono1")

    if spb:
        plist.append("Spb")
        if visu:
            plist.append("Cen")
    if nucleole:
        plist += ['Nuc', 'A_Nuc', 'P_Nuc']

    if telomere:
        plist += ["Telo"]

    snapshot.particles.types = plist

    offset_bond = 0
    offset_particle = 0
    lPolymers = []

    ################################################
    # Polymer chains
    Cen_pos = []
    list_ori = []
    for i in range(Np):

        found_cen = False
        npp = len_chrom[i]  # Number of particles
        # Position of origin of replication
        pos_origins = p_origins[i]

        istrength = None
        if strengths is not None:
            istrength = strengths[i]

        if Sim == []:
            initp = 2 * np.random.rand(3) - 1
        else:
            # print(i)
            initp = Sim.molecules[i].coords[0]

        for p in range(npp - 1):
            inuc = 0
            if nucleole:
                if p in list_nuc[i]:
                    inuc += 1
                if p + 1 in list_nuc[i]:
                    inuc += 1

            snapshot.bonds.group[
                offset_bond + p] = [offset_particle + p,
                                    offset_particle + p + 1]
            if inuc == 0:
                snapshot.bonds.typeid[offset_bond +
                                      p] = bond_list.index('Mono_Mono')  # polymer_A
            if inuc == 1:
                snapshot.bonds.typeid[offset_bond +
                                      p] = bond_list.index('Mono_Nuc')  # polymer_A
            if inuc == 2:
                snapshot.bonds.typeid[offset_bond +
                                      p] = bond_list.index('Nuc_Nuc')  # polymer_A

        offset_bond += npp - 1

        for p in range(npp):
            # print(offset_bond, offset_bond + p)
            if Sim == []:
                new = 2 * (2 * np.random.rand(3) - 1)
                while linalg.norm(initp + new) > R - 1:
                    new = 2 * (2 * np.random.rand(3) - 1)

                initp += new
            else:
                initp = Sim.molecules[i].coords[p]

            snapshot.particles.position[offset_particle + p] = initp

            if p in pos_origins:
                list_ori.append(offset_particle + p)

            if visu and (p in pos_origins):
                snapshot.particles.typeid[
                    offset_particle + p] = plist.index('Ori')  # Ori
            else:
                snapshot.particles.typeid[
                    offset_particle + p] = plist.index('Mono')  # A
                if two_types and p in p_second[i]:
                    snapshot.particles.typeid[
                        offset_particle + p] = plist.index('Mono1')  # A

            if spb and p == Cent[i]:
                Cen_pos.append(offset_particle + p)
                if visu:
                    snapshot.particles.typeid[
                        offset_particle + p] = plist.index('Cen')  # A
                found_cen = True

            if nucleole and p in list_nuc[i]:
                snapshot.particles.typeid[
                    offset_particle + p] = plist.index('Nuc')

            if telomere and (p == 0 or p == npp - 1):
                snapshot.particles.typeid[
                    offset_particle + p] = plist.index('Telo')

        offset_particle += npp

        assert(found_cen == spb)

    phic = 0
    if dump_hic:
        phic = 0 + offset_particle - 1
    ###################################################
    # SPD
    tag_spb = None
    if spb:
        tag_spb = 0 + offset_particle
        # print(tag_spb)
        # print(snapshot.particles[offset_particle])
        snapshot.particles.position[offset_particle] = [-R + 0.1, 0, 0]
        snapshot.particles.typeid[offset_particle] = plist.index('Spb')
        offset_particle += 1

        if not visu:
            for i in range(Np):
                # print(offset_particle - 1, Cen_pos[i])
                snapshot.bonds.group[offset_bond] = [
                    offset_particle - 1, Cen_pos[i]]
                snapshot.bonds.typeid[offset_bond] = bond_list.index(
                    'Spb_Cen')  # polymer_A

                offset_bond += 1
    ############################################################
    # Diffusing elements
    # Defining useful classes

    # Defining particles and bonds for the simulation
    p_tag_list = []
    for i in range(N_diffu):
        npp = 1  # Number of particles

        initp = (R - 2) * (2 * np.random.rand(3) - 1)
        while linalg.norm(initp) > R - 1:
            initp = (R - 2) * (2 * np.random.rand(3) - 1)

        for p in range(npp):

            initp = (R - 2) * (2 * np.random.rand(3) - 1)
            while linalg.norm(initp) > R - 1:
                initp = (R - 2) * (2 * np.random.rand(3) - 1)

            snapshot.particles.position[offset_particle + p] = initp
            snapshot.particles.typeid[
                offset_particle +
                p] = plist.index("Diff")  # Diffu
            p_tag_list.append(offset_particle + p)
        offset_particle += npp

    # Load the configuration

    for i, p in enumerate(snapshot.bonds.group):
        if p[0] == p[1]:
            print(i, p)

    return snapshot, phic, tag_spb, bond_list, plist, Cen_pos, lPolymers, list_ori, p_tag_list


def force_field(traj, bond_list, plist, tag_spb, two_types):

    R = traj["R"]
    micron = traj["micron"]

    # Diffusing elements

    # Yeast case
    spb = traj["spb"]
    nucleole = traj["nucleole"]
    telomere = traj["telomere"]
    microtubule_length = traj["microtubule_length"] * micron
    diameter_nuc = traj["diameter_nuc"] * micron

    # Simulation parameters

    soft = traj["soft"]
    gauss = traj["gauss"]
    assert(type(soft) == bool)
    assert(type(gauss) == bool)

    harmonic = md.bond.harmonic()
    harmonic.bond_coeff.set(bond_list, k=20.0, r0=1)

    harmonic.bond_coeff.set('Mono_Diff', k=10.0, r0=1)

    if spb:
        harmonic.bond_coeff.set('Spb_Cen', k=1000.0, r0=microtubule_length)

    if nucleole:
        harmonic.bond_coeff.set('Nuc_Nuc', k=330, r0=diameter_nuc)
        harmonic.bond_coeff.set(
            'Mono_Nuc', k=330, r0=diameter_nuc / 2. + 1. / 2)

    # Potential for warmup
    if soft:
        def cos_soft(r, rmin, rmax, epsilon, sigma):

            V = epsilon * (1 + np.cos(r * 3.1415 / (rmax)))
            F = epsilon * 3.1415 / (rmax) * np.sin(r * 3.1415 / (rmax))

            return (V, F)

        # nl = md.nlist.tree(r_buff=0.4, check_period=1)
        nl = md.nlist.cell()
        # nl = md.nlist.stencil(r_buff=0.4, check_period=1)
        # nl = md.nlist.cell(r_buff=0.4, check_period=1)

        r_cut = 1.5
        epsilon = 6.5

        table = md.pair.table(width=1000, nlist=nl)
        table.pair_coeff.set(plist, plist,
                             func=cos_soft, rmin=0, rmax=r_cut,
                             coeff=dict(epsilon=epsilon, sigma=1.0))

        if nucleole:
            for ip1, p1 in enumerate(plist):
                for p2 in plist[ip1:]:
                    inuc = 0
                    if "Nuc" in p1:
                        inuc += 1
                    if "Nuc" in p2:
                        inuc += 1
                    if inuc == 1:
                        d = 0.5 + diameter_nuc / 2.
                        d = r_cut * d
                    if inuc == 2:
                        d = r_cut * diameter_nuc
                        # smaller here
                    if inuc == 0:
                        continue
                    table.pair_coeff.set(p1, p2,
                                         func=cos_soft, rmin=0, rmax=d,
                                         coeff=dict(epsilon=epsilon, sigma=d))

    else:

        if gauss:
            r_cut = 1.5
            # nl = md.nlist.tree(r_buff=0.4, check_period=1)
            nl = md.nlist.cell()

            # gauss = md.pair.gauss(r_cut=r_cut, nlist=nl)
            # gauss.pair_coeff.set(plist, plist, epsilon=1.0, sigma=0.3)

            def gauss_center_decay_strength(r, rmin, rmax, c=0, sigma=0.3, epsilon=1):

                V = epsilon * np.exp(-(r - c)**2 / (2 * sigma**2))
                F = epsilon * (r - c) / sigma**2 * np.exp(-(r - c)**2 / (2 * sigma**2))

                return (V, F)

            table = md.pair.table(width=1000, nlist=nl)
            table.pair_coeff.set(plist, plist,
                                 func=gauss_center_decay_strength, rmin=0, rmax=1.5,
                                 coeff=dict(epsilon=1, sigma=.3))
            if two_types:
                def gauss_center_decay_strength_a(r, rmin, rmax, c=0, sigma=0.3, epsilon=1, epsilona=traj.get("epsilona", -0.2)):

                    V1, F1 = gauss_center_decay_strength(
                        r, rmin, rmax, c=c, sigma=sigma, epsilon=epsilon)
                    Va, Fa = gauss_center_decay_strength(
                        r, rmin, rmax, c=sigma * 1.6, sigma=sigma / 2, epsilon=epsilona)
                    return (V1 + Va, F1 + Fa)
                table.pair_coeff.set(["Mono1"], ["Mono1"],
                                     func=gauss_center_decay_strength_a, rmin=0, rmax=1.5,
                                     coeff=dict(epsilon=1, sigma=.3))

        else:
            r_cut = 1.12
            # nl = md.nlist.tree()  # r_buff=10, check_period=1)
            nl = md.nlist.cell()
            gauss = md.pair.lj(r_cut=r_cut, nlist=nl)  # , d_max=diameter_nuc)

            gauss.pair_coeff.set(plist, plist, epsilon=1.0, sigma=0.3)

        if nucleole:
            for ip1, p1 in enumerate(plist):
                for p2 in plist[ip1:]:
                    inuc = 0
                    if "Nuc" in p1:
                        inuc += 1
                    if "Nuc" in p2:
                        inuc += 1
                    if inuc == 1:
                        gauss.pair_coeff.set(
                            p1,
                            p2,
                            epsilon=.5,
                            sigma=0.5 +
                            diameter_nuc /
                            2.,
                            r_cut=(
                                0.5 +
                                diameter_nuc /
                                2.) * r_cut)
                    if inuc == 2:
                        gauss.pair_coeff.set(p1, p2, epsilon=1.0, sigma=diameter_nuc,
                                             r_cut=diameter_nuc * r_cut)
    # gauss.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
    # gauss.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)

    # Spherical confinement
    sphere = md.wall.group()
    r_extrap = 0.95
    sphere.add_sphere(r=R, origin=(0.0, 0.0, 0.0), inside=True)
    # lj much more slower (at least in thu minimisation)
    wall_force_slj = md.wall.lj(sphere, r_cut=1.12)
    wall_force_slj.force_coeff.set(plist, epsilon=1.0, sigma=1.0,
                                   r_cut=1.12, mode="shift", r_extrap=r_extrap)

    if spb:
        wall_force_slj.force_coeff.set("Spb", epsilon=1.0, sigma=1.0,
                                       r_cut=-1, mode="shift")

    # wall_force_slj.set_params(mode="shift")

    if nucleole:
        wall_force_slj.force_coeff.set(
            'Nuc',
            epsilon=1.0,
            sigma=diameter_nuc,
            r_cut=diameter_nuc * 1.12, mode="shift", r_extrap=diameter_nuc * r_extrap)
    if telomere:
        wall_force_slj.force_coeff.set("Telo", epsilon=2.0, sigma=1.5,
                                       r_cut=3, mode="shift", r_extrap=r_extrap)

    # Group;
    all_beads = group.all()
    Spb_g = None
    if spb:
        Spb_g = group.tag_list(name="Spb", tags=[tag_spb])
        pspb = [p.position for p in Spb_g]
        print(pspb)

        all_move = group.difference(name="move", a=all_beads, b=Spb_g)
    else:
        all_move = all_beads

    return all_beads, all_move, Spb_g, nl


def minimize(traj, all_move, system, snapshot, Spb_g, Cen_pos, microtubule_length):

    R = traj["R"]
    data_folder = traj["data_folder"]

    # Diffusing elements

    # Yeast case
    spb = traj["spb"]

    visu = traj["visu"]

    # Scenari

    converged = False
    dt = 0.002

    dcd = dump.dcd(filename=data_folder + 'init.dcd',
                   period=10, overwrite=True)
    while not converged and not visu:
        try:
            if '2.2' in hoomd.__version__:
                method = md.integrate.mode_minimize_fire(dt=dt)
                nve = md.integrate.nve(group=all_move)
            else:
                method = md.integrate.mode_minimize_fire(group=all_move, dt=dt)
            while not(method.has_converged()):

                if spb:
                    pspb = [p.position for p in Spb_g]

                    for cen in Cen_pos:
                        cent_tmp = system.particles[cen]
                        d = linalg.norm(
                            np.array(pspb[0]) - np.array(cent_tmp.position))
                        if d > 2 * microtubule_length:
                            print("MT too long", d)
                            raise

                    """
                    print(pspb)
                    for cen in Cen_pos:
                        cent_tmp = system.particles[cen]
                        print(cent_tmp.position)
                        print(linalg.norm(np.array(pspb[0])-np.array(cent_tmp.position)))
                        print(R * microtubule_length)
                    """
                # exit()
                hoomd.run(100)

                # Check all particles are in the sphere:
                for p in system.particles:
                    if linalg.norm(p.position) > R:
                        raise

            converged = True
            if '2.2' in hoomd.__version__:
                nve.disable()

        except:
            if '2.2' in hoomd.__version__:
                nve.disable()
            converged = False
            dt /= 2.
            print("Reducing time step", dt)
            # Restore positions
            for ip, p in enumerate(snapshot.particles.position):

                system.particles[ip].position = p

    dcd.disable()


def simulate(traj):

    seed = traj["seed"]

    micron = traj["micron"]
    data_folder = traj["data_folder"]

    # Diffusing elements
    cut_off_inte = traj["cut_off_inte"]
    p_inte = traj["p_inte"]
    p_off = traj["p_off"]
    sim_dt = traj["sim_dt"]

    dscale = traj["dscale"]
    # Yeast case
    spb = traj["spb"]
    assert(type(spb) == bool)

    microtubule_length = traj["microtubule_length"] * micron
    visu = traj["visu"]
    dump_hic = traj["dump_hic"]
    two_types = traj.get("two_types", False)

# Scenari
    diff_alone = traj["diff_alone"]
    diff_bind_when_free = traj["diff_bind_when_free"]
    diff_bind_when_on_DNA = traj["diff_bind_when_on_DNA"]
    assert(type(diff_alone) == bool)
    assert(type(diff_bind_when_on_DNA) == bool)
    assert(type(diff_bind_when_free) == bool)

    # Simulation parameters
    n_steps = traj["n_steps"]
    length_steps = traj["length_steps"]
    benchmark = traj["benchmark"]
    warmup = traj["warmup"]
    dcd_period = traj["dcd_period"]

    np.random.seed(seed)
    hoomd.context.initialize()  # "--mode=cpu ")

    if diff_alone:
        # Check
        assert(diff_bind_when_free is False)
        assert (diff_bind_when_on_DNA is False)

    # End of parameter
    ##########################################

    #########################################
    # Define polymer bonding and positions

    snapshot, phic, tag_spb, bond_list, plist, Cen_pos, lPolymers, list_ori, p_tag_list = \
        create_initial_configuration(traj)
    system = init.read_snapshot(snapshot)

    for i, p in enumerate(system.particles):
        # print(p)
        # exit()
        assert p.tag == i

    for i, b in enumerate(system.bonds):
        if b.a == b.b:
            print(b.a, b.b)

            raise
        assert b.tag == i
    ###############################################

    ###############################################
    # Defining force field:
    all_beads, all_move, Spb_g, nl = force_field(
        traj, bond_list=bond_list, plist=plist, tag_spb=tag_spb, two_types=two_types)

    # Log
    if not visu:
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

    # Warmup

    minimize(traj, all_move, system, snapshot, Spb_g, Cen_pos, microtubule_length)

    # Dumping

    if visu:
        xml = deprecated.dump.xml(
            filename=data_folder +
            "atoms.hoomdxml",
            period=None,
            group=all_beads,
            vis=True)
        # xml.disable()
        return
    # gsd = dump.gsd(filename=data_folder + "atoms.gsd",period=None,group=all_beads)

    # Dynamics

    def Change_type(typep, particle_list, snp, msg=""):
        # print(particle_list)
        for p in particle_list:
            if "Ori" in typep:
                # Remove it from the list activated
                # list_ori.remove(p)
                pass
            else:
                snp.particles[p].type = typep
        if particle_list != [] and msg != "":
            print(msg)

    def Bind(typeb, bondlist, snp):
        btags = []
        for b1, b2 in bondlist:
            btags.append(snp.bonds.add(typeb, b1, b2))
        return btags

    def Release(btags, snp):
        for bt in btags:
            snp.bonds.remove(bt)

    def AddParticle(position, type):
        snp.particles.add(type)
        snp.particles[-1].position = position

    def Shift(bonds, snp):
        for tag, new in bonds:
            b = snp.bonds.get(tag)
            btype = "" + b.type
            fork = b.b + 0
            snp.bonds.remove(tag)

            # print(b.type)
            snp.bonds.add(btype, new, fork)
            # print(new,b)
            # print(dir(snp.bonds))
            # b.a = new

    snp = system  # .take_snapshot()

    # nl.tune(warmup=1,steps=1000)

    # Small warmup

    t0 = time.time()
    md.integrate.mode_standard(dt=sim_dt)
    method = md.integrate.langevin(group=all_move, kT=1, seed=seed, dscale=dscale)

    if benchmark:
        print(nl.tune(warmup=4000, r_min=0.3, r_max=0.8, jumps=5, steps=5000))
        return

    md.integrate.mode_standard(dt=sim_dt / 4)
    hoomd.run(100)
    md.integrate.mode_standard(dt=sim_dt / 2)
    hoomd.run(100)
    md.integrate.mode_standard(dt=sim_dt)

    if warmup != 0:
        hoomd.run(warmup)

    dcd = dump.dcd(filename=data_folder + 'poly.dcd',
                   period=dcd_period, overwrite=True)

    r_hic = []
    if dump_hic:
        group_hic = group.tags(name="hic", tag_min=0, tag_max=phic)

    global timeit
    timeit = True
    global t0
    t0 = time.time()

    def Timeit(where=""):
        global timeit
        global t0
        if timeit:
            if where == "":
                print(time.time() - t0)
            else:
                print(where, time.time() - t0)
            t0 = time.time()

    bonds = []
    for i in range(n_steps):

        # Chek that the microtubule length is correct
        if spb:
            for cen in Cen_pos:
                cent_tmp = system.particles[cen]
                # print(cent_tmp.position)
                pspb = [p.position for p in Spb_g]
                d = linalg.norm(
                    np.array(pspb[0]) - np.array(cent_tmp.position))
                if d > 2 * microtubule_length:
                    print("MT too long", d)
                    exit()

        # Dump the Hi-Cs
        Timeit()
        # system.restore_snapshot(snp)

        hoomd.run(length_steps // 2, profile=False)
        Timeit("After first half")

        if dump_hic and i % traj.get("hic_period", 1) == 0:
            ph = np.array([p.position for p in group_hic])

            D = cdist(ph, ph)
            D[D < 2] = 1
            D[D >= 2] = 0
            np.fill_diagonal(D, 0)
            if r_hic != []:
                r_hic += D
            else:
                r_hic = D
            if i % 32 == 0:
                np.save(data_folder + "/hic", r_hic)

        # snp = system.take_snapshot()

        # Bond tags to release (Alone particle)
        remv = []
        for ib, b in enumerate(bonds):
            if np.random.rand() > p_off:
                continue

            p_tag_list.append(b[2])
            list_ori.append(b[1])

            Release([b[0]], snp)
            remv.append(ib)

        remv.sort()
        for ib in remv[::-1]:
            bonds.pop(ib)

        Timeit("AFter update")
        hoomd.run(length_steps // 2, profile=True)
        Timeit("AFter second half")

        # Update Type because of (Ori to passivated)

        # Update group

        # Find new interacting particles

        # First check if Dimer are close from one origin

        print("LAAAAAAAAAAAAAAA", p_tag_list, list_ori, bonds)
        p_diffu = np.array([snp.particles[diff].position for diff in p_tag_list])

        p_origin = np.array([snp.particles[ori].position for ori in list_ori])

        if not(p_tag_list != [] and list_ori != []):
            print("No interactions")

        if p_tag_list != [] and list_ori != []:
            distances = cdist(p_diffu, p_origin)
            print(distances.shape)
            # Reorder the distances with the dimer tags

            # Groups Diff-Diff by bond to compute the distances

            activated_ori = []
            activated_p = []

            for iD, ptag in enumerate(p_tag_list):
                # print(d2.shape)
                # print(d2[iD])
                lo = list(range(len(distances[iD])))
                np.random.shuffle(lo)

                for iorigin in lo:
                    di = distances[iD][iorigin]
                    if iorigin in activated_ori:
                        # Needed because we don't want an origin to be activated
                        # twice
                        continue
                    if di > cut_off_inte:
                        continue

                    if np.random.rand() > p_inte:
                        continue
                    activated_ori.append(iorigin + 0)
                    activated_p.append(iD)

                    new_btags = Bind("Mono_Diff", [[list_ori[iorigin], ptag]], snp)

                    bonds.append([new_btags[0], list_ori[iorigin], ptag])

                    break
            activated_ori.sort()
            for ori in activated_ori[::-1]:
                list_ori.pop(ori)

            activated_p.sort()
            for p in activated_p[::-1]:
                p_tag_list.pop(p)

        Timeit("After binding")
        # t0 = time.time()

        # print(time.time() -t0)
        # Then if it is the case attach them according to p law to the origin

    # print(gauss.get_energy(all_beads), wall_force_slj.get_energy(all_beads))
    print(time.time() - t0)
    logger.disable()
    method.disable()
    dcd.disable()


if __name__ == "__main__":
    import sys

    traj = load_parameters(sys.argv[1])
    simulate(traj)
