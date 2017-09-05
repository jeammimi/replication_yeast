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
    if p_second != []:
        ps = []
        for ch in p_second:
            ps.append([])
            for p1, p2 in ch:
                ps[-1] += range(p1, p2)
        p_second = ps

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

        lPolymers.append(Polymer(i,
                                 offset_particle,
                                 offset_particle + npp - 1,
                                 [po + offset_particle for po in pos_origins]))
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
        npp = 2  # Number of particles

        initp = (R - 2) * (2 * np.random.rand(3) - 1)
        while linalg.norm(initp) > R - 1:
            initp = (R - 2) * (2 * np.random.rand(3) - 1)
        if diff_bind_when_free:
            for p in range(npp - 1):
                snapshot.bonds.group[
                    offset_bond +
                    p] = [
                    offset_particle +
                    p,
                    offset_particle +
                    p +
                    1]
                snapshot.bonds.typeid[offset_bond +
                                      p] = bond_list.index('Diff_Diff')  # Diff_Diff
            offset_bond += npp - 1

        p_tag_list.append([])
        for p in range(npp):

            # print(offset_bond, offset_bond + p)
            if diff_bind_when_free:
                new = 2 * (2 * np.random.rand(3) - 1)
                while linalg.norm(initp + new) > R - 1:
                    new = 2 * (2 * np.random.rand(3) - 1)
                    # print(initp,new,R,linalg.norm(initp + new))
                    # exit()
                initp += new
            else:
                initp = (R - 2) * (2 * np.random.rand(3) - 1)
                while linalg.norm(initp) > R - 1:
                    initp = (R - 2) * (2 * np.random.rand(3) - 1)

            snapshot.particles.position[offset_particle + p] = initp
            snapshot.particles.typeid[
                offset_particle +
                p] = plist.index("Diff")  # Diffu
            p_tag_list[-1].append(offset_particle +
                                  p)
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
        except:
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
    fork_speed = traj["fork_speed"]
    dt_speed = traj["dt_speed"]
    dscale = traj["dscale"]
    # Yeast case
    spb = traj["spb"]
    assert(type(spb) == bool)

    microtubule_length = traj["microtubule_length"] * micron
    Activ_Origins = traj["Activ_Origins"]
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

    ramp_type = traj.get("ramp_type", "exp")
    ramp = traj.get("ramp_time", 3)

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

    ramp_type = "exp"
    ramp = 3
    if ramp_type == "exp":
        for couple in p_tag_list:
            Change_type("I_Diff", couple, snp)

    group_diffu = group.type(name="Diff", type='Diff')

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
    Ndiff_libre_t = []

    N_diffu = traj["N_diffu"]

    offset_diff = np.min(p_tag_list)
    print(offset_diff)

    record_diffusing = [Diffusing(d) for d in np.arange(N_diffu * 2)]

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

    previous_actifs = 0
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

        N_actifs = int(len(p_tag_list) * (1 - np.exp(- i * dt_speed / ramp)))
        print(previous_actifs, N_actifs)
        for couple in p_tag_list[previous_actifs:N_actifs]:
            Change_type("Diff", couple, snp)
            print("Activated", couple)
        previous_actifs = N_actifs

        hoomd.run(length_steps // 2, profile=False)
        Timeit("After first half")

        if dump_hic:
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

        # update the position of the monomer by updating bonds
        ended = 0
        for iP, P in enumerate(lPolymers):
            verbose = False
            # if iP == 9:
            #    verbose = True
            bind_diff, diff_diff, shifted_bonds, \
                passivated_origin, to_release, alone = P.increment_time(
                    dt_speed, verbose)

            ###################################################################
            # Only to keep track of the diffusing elements
            for diff1, diff2 in bind_diff:

                record_diffusing[
                    diff1 - offset_diff].end_replication(i * dt_speed, pos=snp.particles[diff1].position)
                record_diffusing[
                    diff2 - offset_diff].end_replication(i * dt_speed, pos=snp.particles[diff2].position)

            for diff in alone:
                if record_diffusing[diff - offset_diff].bound:
                    record_diffusing[
                        diff - offset_diff].end_bound(i * dt_speed, pos=snp.particles[diff].position)
                elif record_diffusing[diff - offset_diff].replicating:
                    record_diffusing[
                        diff - offset_diff].end_replication(i * dt_speed, pos=snp.particles[diff].position)
                else:
                    print(diff, record_diffusing[diff - offset_diff].free)
                    raise

            ###################################################################

            ###################################################################
            # take care of the bondings

            for ori in passivated_origin:
                list_ori.remove(ori)

            if not diff_alone:
                Shift(shifted_bonds, snp)
                # Bond tags to release (Alone particle)
                Release(to_release, snp)

                if diff_bind_when_free:
                    # Pair of diffu to attach
                    Bind("Diff_Diff", bind_diff, snp)
                    # We cannot use the single diff anymore
                    Change_type("S_Diff", alone, snp)
                    # Change type for pair of diff diff
                    Change_type("Diff", diff_diff, snp)

            # check for releasing alone binded elements

            for ori_not_started in P.get_free_origins():
                diff = P.get_diff_at_origin(ori_not_started)
                if diff != []:
                    if np.random.rand() > p_off:
                        continue
                    ptag, bond_tag = diff[0]
                    P.dettach_one_diff(ptag, ori_not_started)
                    Release([bond_tag], snp)
                    Change_type("Diff", [ptag], snp)

            if P.modules == []:
                ended += 1

        Timeit("AFter update")
        hoomd.run(length_steps // 2, profile=True)
        Timeit("AFter second half")

        group_diffu.force_update()
        # Update Type because of (Ori to passivated)

        # Update group

        # Find new interacting particles

        # First check if Dimer are close from one origin

        p_diffu = np.array([p.position for p in group_diffu])
        tag_diffu = [p.tag for p in group_diffu]

        p_origin = np.array([snp.particles[ori].position for ori in list_ori])

        Ndiff_libre_t.append(len(tag_diffu))

        if not(tag_diffu != [] and list_ori != []):
            print("No interactions")

        if tag_diffu != [] and list_ori != []:
            distances = cdist(p_diffu, p_origin)
            print(distances.shape)
            # Reorder the distances with the dimer tags
            Indexes = []
            PTags = []
            # t0 = time.time()
            Btags = []
            # Groups Diff-Diff by bond to compute the distances

            if diff_bind_when_free:
                for b in system.bonds:
                    if b.type == 'Diff_Diff' and system.particles[
                            b.a].type == 'Diff':
                        if b.a in tag_diffu:
                            Indexes.append(tag_diffu.index(b.a))
                            Indexes.append(tag_diffu.index(b.b))
                            Btags.append(b.tag)
                            PTags.append([b.a, b.b])

                # print(time.time() -t0)

                d2 = distances[Indexes][::2] / 2 + distances[Indexes][1::2] / 2
            else:
                n_diffu = len(tag_diffu)
                Indexes = list(range(n_diffu))
                Btags = [None] * n_diffu
                PTags = [[t] for t in tag_diffu]
                d2 = distances[Indexes]

            activated = []
            for iD, (btag, ptags) in enumerate(zip(Btags, PTags)):
                # print(d2.shape)
                # print(d2[iD])
                for iorigin, di in enumerate(d2[iD]):
                    if iorigin in activated:
                        # Needed because we don't want an origin to be activated
                        # twice
                        continue
                    if di > cut_off_inte:
                        continue
                    if np.random.rand() > p_inte:
                        continue

                    for P in lPolymers:
                        if not P.has_origin(list_ori[iorigin]):
                            continue

                        if diff_bind_when_free and \
                           not diff_bind_when_on_DNA:
                            Release([btag], snp)  # Break the dimer
                            btag = None  # We need btag only in the case where they stays attached

                        if not diff_alone:
                            # Or attached separatly or already bound:

                            if diff_bind_when_free:

                                # We are sure they are two and we can
                                # start
                                Change_type(
                                    'F_Diff', ptags, snp)  # Diffusive element attached
                                particular_origin = list_ori[iorigin]
                                new_btags = Bind("Mono_Diff", [[particular_origin, ptags[0]], [
                                                 particular_origin, ptags[1]]], snp)

                                activated.append(0 + iorigin)
                                P.add_fork(
                                    ptags, particular_origin, new_btags, btag, fork_speed=fork_speed)

                                for diff in ptags:
                                    record_diffusing[
                                        diff - offset_diff].start_replication(particular_origin, i * dt_speed, pos=snp.particles[diff].position)

                            else:
                                Change_type(
                                    'F_Diff', ptags, snp)  # Diffusive element attached
                                particular_origin = list_ori[iorigin]
                                new_btags = Bind(
                                    "Mono_Diff", [[particular_origin, ptags[0]]], snp)
                                start = P.attach_one_diff(
                                    ptags[0], particular_origin, new_btags[0])

                                if start:
                                    # get particles involves
                                    p1, p2 = P.get_diff_at_origin(
                                        particular_origin)
                                    if diff_bind_when_on_DNA:
                                        btag = Bind("Diff_Diff", [
                                            [p1[0], p2[0]]], snp)[0]

                                    activated.append(0 + iorigin)
                                    P.add_fork([p1[0], p2[0]], particular_origin,
                                               [p1[1], p2[1]], btag, fork_speed=fork_speed)

                                    record_diffusing[
                                        p1[0] - offset_diff].start_replication(particular_origin, i * dt_speed, pos=snp.particles[p1[0]].position)
                                    record_diffusing[
                                        p2[0] - offset_diff].start_replication(particular_origin, i * dt_speed, pos=snp.particles[p2[0]].position)
                                else:
                                    record_diffusing[ptags[0] -
                                                     offset_diff].start_bound(particular_origin, i * dt_speed, pos=snp.particles[ptags[0]].position)

                        else:
                            # start when touched and release
                            particular_origin = list_ori[iorigin]
                            activated.append(iorigin)

                            P.add_fork([None, None], particular_origin, [
                                       None, None], None, fork_speed=fork_speed)

                        break
                    # If we arrive there it means that one interaction has beeen
                    # found
                    break
            activated.sort()
            print(activated)
            print(list_ori)

            for io in activated[::-1]:
                print(io)
                list_ori.pop(io)
        Timeit("After binding")
        # t0 = time.time()
        with open(data_folder + "polymer_timing.dat", "wb") as f:
            cPickle.dump(lPolymers, f, protocol=2)
        with open(data_folder + "Ndiff_libre_t.dat", "wb") as f:
            cPickle.dump(Ndiff_libre_t, f, protocol=2)

        with open(data_folder + "record_diffusing.dat", "wb") as f:
            cPickle.dump(record_diffusing, f, protocol=2)

        Timeit("After writing")

        if traj.get("early_stop", False) and list_ori == [] and ended == len(lPolymers):
            break

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
