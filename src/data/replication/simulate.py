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
from PMotion import Polymer
import _pickle as cPickle
from createPoly import create_init_conf_yeast
import time
import json


def simulate(traj_filename):

    with open(traj_filename, "r") as f:
        traj = json.load(f)

    seed = traj["seed"]
    len_chrom = traj["len_chrom"]
    Cent = traj["Cent"]
    p_ribo = traj["p_ribo"]
    R = traj["R"]
    micron = traj["micron"]
    data_folder = traj["data_folder"]

    # Diffusing elements
    N_diffu = traj["N_diffu"]
    cut_off_inte = traj["cut_off_inte"]
    p_inte = traj["p_inte"]
    dt = traj["dt"]
    p_origins = traj["p_origins"]

    # Yeast case
    spb = traj["spb"]
    nucleole = traj["nucleole"]
    telomere = traj["telomere"]
    microtubule_length = traj["microtubule_length"] * micron
    diameter_nuc = traj["diameter_nuc"] * micron
    special_start = traj["special_start"]
    Activ_Origins = traj["Activ_Origins"]
    visu = traj["visu"]
    dump_hic = traj["dump_hic"]

    # Scenari
    diff_alone = traj["diff_alone"]
    diff_bind_when_free = traj["diff_bind_when_free"]
    diff_bind_when_on_DNA = traj["diff_bind_when_on_DNA"]
    replicate_DNA = traj["replicate_DNA"]


    np.random.seed(seed)
    hoomd.context.initialize("--mode=cpu")


    if diff_alone:
        # Check
        assert(diff_bind_when_free is False)
        assert (diff_bind_when_on_DNA is False)

    # End of parameter
    ##########################################

    #########################################
    # Define polymer bonding and positions

    Np = len(len_chrom)
    assert(len(len_chrom) == len(Cent) == len(p_ribo))
    if special_start:
        Sim = create_init_conf_yeast(
            len_chrom=len_chrom,
            dist_centro=Cent,
            p_ribo=p_ribo,
            Radius=R,
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

    plist = ['Mono', 'Ori', 'Diff', 'A_Ori', 'P_Ori', 'S_Diff', 'F_Diff']

    if spb:
        plist.append("Spb")
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
                offset_bond +
                p] = [
                offset_particle +
                p,
                offset_particle +
                p +
                1]
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
                snapshot.particles.typeid[
                    offset_particle + p] = plist.index('Ori')  # Ori
            else:
                snapshot.particles.typeid[
                    offset_particle + p] = plist.index('Mono')  # A

            if spb and p == Cent[i]:
                Cen_pos.append(offset_particle + p)

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
                initp = (R - 1) * (2 * np.random.rand(3) - 1)

            snapshot.particles.position[offset_particle + p] = initp
            snapshot.particles.typeid[
                offset_particle +
                p] = plist.index("Diff")  # Diffu

        offset_particle += npp

    # Load the configuration

    for i, p in enumerate(snapshot.bonds.group):
        if p[0] == p[1]:
            print(i, p)

    system = init.read_snapshot(snapshot)

    for i, p in enumerate(system.particles):
        # print(p)
        # exit()
        assert p.tag == i

    for i, b in enumerate(system.bonds):
        if b.a == b.b:
            print(b.a, b.b)

            raise
        # print(p)
        # exit()
        assert b.tag == i
    ###############################################

    ###############################################
    # Defining force field:
    harmonic = md.bond.harmonic()
    harmonic.bond_coeff.set(bond_list, k=330.0, r0=1)

    harmonic.bond_coeff.set('Mono_Diff', k=10.0, r0=1)

    if spb:
        harmonic.bond_coeff.set('Spb_Cen', k=1000.0, r0=microtubule_length)

    if nucleole:
        harmonic.bond_coeff.set('Nuc_Nuc', k=330, r0=diameter_nuc)
        harmonic.bond_coeff.set(
            'Mono_Nuc', k=330, r0=diameter_nuc / 2. + 1. / 2)

    nl = md.nlist.tree(r_buff=0.4, check_period=1)

    # Potential for warmup
    gauss = md.pair.gauss(r_cut=3.0, nlist=nl)

    gauss.pair_coeff.set(plist, plist, epsilon=1.0, sigma=1.0)

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
                            2.) *
                        3)
                if inuc == 2:
                    gauss.pair_coeff.set(p1, p2, epsilon=1.0, sigma=diameter_nuc,
                                         r_cut=3 * diameter_nuc)
    # gauss.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
    # gauss.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)

    # Spherical confinement
    sphere = md.wall.group()
    sphere.add_sphere(r=R, origin=(0.0, 0.0, 0.0), inside=True)

    wall_force_slj = md.wall.slj(sphere, r_cut=3.0)
    wall_force_slj.force_coeff.set(plist, epsilon=1.0, sigma=1.0, r_cut=1.12)

    if nucleole:
        wall_force_slj.force_coeff.set(
            'Nuc',
            epsilon=1.0,
            sigma=diameter_nuc,
            r_cut=diameter_nuc * 1.12)
    if telomere:
        wall_force_slj.force_coeff.set(plist, epsilon=2.0, sigma=1.5, r_cut=3)

    # Group;
    all_beads = group.all()
    if spb:
        Spb_g = group.tag_list(name="Spb", tags=[tag_spb])
        pspb = [p.position for p in Spb_g]
        print(pspb)

        all_move = group.difference(name="move", a=all_beads, b=Spb_g)
    else:
        all_move = all_beads
    # Log
    logger = analyze.log(
        filename=data_folder +
        'mylog.log',
        period=1000,
        quantities=[
            'temperature',
            'potential_energy',
            'kinetic_energy',
            'volume',
            'pressure'],
        overwrite=True)

    # Warmup
    converged = False
    dt = 0.005
    while not converged and not visu:
        try:

            method = md.integrate.mode_minimize_fire(group=all_move, dt=dt)
            while not(method.has_converged()):

                if spb:
                    pspb = [p.position for p in Spb_g]
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
            converged = True
        except:
            converged = False
            dt /= 2.
            print(dt)
            # Restore positions
            for ip, p in enumerate(snapshot.particles.position):

                system.particles[ip].position = p

    """
    gauss.disable()

    slj=md.pair.slj(r_cut=2, nlist=nl)
    slj.pair_coeff.set(plist,plist,sigma=1,epsilon=1,r_cut=1.12)
    print("Second minimizing")
    method=md.integrate.mode_minimize_fire(group=all_beads,dt=0.05)
    while not(method.has_converged()):
       hoomd.run(100)
    """
    # hoomd.run(1000000)
    # method.disable()

    # Dumping

    if visu:
        xml = deprecated.dump.xml(
            filename=data_folder +
            "atoms.hoomdxml",
            period=None,
            group=all_beads,
            vis=True)
        exit()
    # gsd = dump.gsd(filename=data_folder + "atoms.gsd",period=None,group=all_beads)
    dcd = dump.dcd(filename=data_folder + 'poly.dcd',
                   period=100, overwrite=True)

    # Dynamics

    t0 = time.time()
    md.integrate.mode_standard(dt=0.01)
    method = md.integrate.langevin(group=all_move, kT=1, seed=seed)
    snp = system  # .take_snapshot()

    def Change_type(typep, particle_list, snp, msg=""):
        # print(particle_list)
        for p in particle_list:
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

    group_diffu = group.type(name="Diff", type='Diff')

    if Activ_Origins != []:
        group_origin = group.type(name="Activ_Ori", type=Activ_Origins[0])
        if len(Activ_Origins) > 1:
            for t in Activ_Origins[1:]:
                group_origin = group.union(name="Activ_origin", a=group_origin,
                                           b=group.type(name="tmp", type=t))

    r_hic = []
    if dump_hic:
        group_hic = group.tags(name="hic", tag_min=0, tag_max=phic)
    # nl.tune(warmup=1,steps=1000)

    for i in range(100):

        # Chek that the microtubule length is correct
        if spb:
            for cen in Cen_pos:
                cent_tmp = system.particles[cen]
                # print(cent_tmp.position)
                d = linalg.norm(
                    np.array(pspb[0]) - np.array(cent_tmp.position))
                if d > 2 * microtubule_length:
                    print("MT too long", d)
                    exit()

        # Dump the Hi-Cs

        # system.restore_snapshot(snp)
        hoomd.run(1000)

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
            np.save(data_folder + "/hic", r_hic)

        # snp = system.take_snapshot()

        # update the position of the monomer by updating bonds

        for iP, P in enumerate(lPolymers):
            verbose = False
            # if iP == 9:
            #    verbose = True
            bind_diff, diff_diff, shifted_bonds, \
                passivated_origin, to_release, alone = P.increment_time(
                    1, verbose)

            Change_type(
                'P_Ori',
                passivated_origin,
                snp,
                msg="")  # Passivated origin

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

        group_diffu.force_update()
        group_origin.force_update()
        # Update Type because of (Ori to passivated)

        # Update group

        # Find new interacting particles

        # First check if Dimer are close from one origin

        p_diffu = np.array([p.position for p in group_diffu])
        tag_diffu = [p.tag for p in group_diffu]

        p_origin = np.array([p.position for p in group_origin])
        tag_origin = [p.tag for p in group_origin]

        if tag_diffu != [] and tag_origin != []:
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
                        if not P.has_origin(tag_origin[iorigin]):
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
                                particular_origin = tag_origin[iorigin]
                                new_btags = Bind("Mono_Diff", [[particular_origin, ptags[0]], [
                                                 particular_origin, ptags[1]]], snp)
                                Change_type(
                                    'A_Ori', [particular_origin], snp)
                                activated.append(iorigin)
                                P.add_fork(
                                    ptags, particular_origin, new_btags, btag)

                            else:
                                Change_type(
                                    'F_Diff', ptags, snp)  # Diffusive element attached
                                particular_origin = tag_origin[iorigin]
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

                                    Change_type(
                                        'A_Ori', [particular_origin], snp)
                                    P.add_fork([p1[0], p2[0]], particular_origin,
                                               [p1[1], p2[1]], btag)

                        else:
                            # start when touched and release
                            particular_origin = tag_origin[iorigin]
                            activated.append(iorigin)
                            Change_type(
                                'A_Ori', [particular_origin], snp)
                            P.add_fork([None, None], particular_origin, [
                                       None, None], None)

                        break
                    # If we arrive there it means that one interaction has beeen
                    # found
                    break
        # t0 = time.time()
        with open(data_folder + "polymer_timing.dat", "wb") as f:
            cPickle.dump(lPolymers, f)
        # print(time.time() -t0)
        # Then if it is the case attach them according to p law to the origin

    print(gauss.get_energy(all_beads), wall_force_slj.get_energy(all_beads))
    print(time.time() - t0)


if __name__ == "__main__":
    simulate("./default.json")
