import pandas
import json


def load_ori_position(File, ori_type, lengths, coarse, verbose=True, strength=None, coarsed=True):

    Oris = pandas.read_csv(File, comment="#")

    l_ori = [[] for i in range(len(lengths))]
    strengths = [[] for i in range(len(lengths))]

    # Filters ori:
    istrength = [[] for i in range(len(lengths))]
    for ch, p, status in zip(Oris["chr"], Oris["start"], Oris["status"]):
        if status in ori_type:  # ,"Likely"]:
            if coarsed:
                l_ori[ch - 1].append(int(p / coarse))
            else:
                l_ori[ch - 1].append(p / coarse)

            if strength is not None:
                strengths[ch - 1].append(strength[status])
            else:
                strengths[ch - 1].append(1)
            istrength[ch - 1].append(status)

    # Then remove duplicate (/kb ) and outside of boundaries

    tot = 0
    for i in range(len(lengths)):
        isize = len(l_ori[i])
        if coarsed:
            l_ori[i] = list(set(l_ori[i]))
        l_ori[i].sort()
        if verbose:
            print(isize, len(l_ori[i]))
        while not (max(l_ori[i]) <= lengths[i] + 0.99):
            print(max(l_ori[i]), lengths[i])
            l_ori[i].remove(max(l_ori[i]))
            strengths[i].pop(-1)
            istrength[i].pop(-1)

        # print(max(l_ori[i]),len_chrom[i]*5

        tot += len(l_ori[i])
        assert(max(l_ori[i]) < lengths[i] + 0.99)

    if strength is not None:
        return l_ori, strengths, istrength
    return l_ori


def load_lengths_and_centro(File, coarse=1, verbose=True):
    gff = pandas.read_csv(File, sep="\t", comment="#", header=None, names=[
                          "chr", "SGD", "info", "start", "end", "m1", "m2", "m3", "comment"], low_memory=False)

    lengths = [int(end) // int(coarse) for inf, end, chro in zip(gff["info"], gff["end"], gff["chr"])
               if inf == "chromosome" and "mt" not in chro]
    # for i in range(remove):
    # print(lengths,len(lengths))
    cents = [int(end) // int(coarse)
             for inf, end in zip(gff["info"], gff["end"]) if inf == "centromere"]
    # print(cents,len(cents))
    if verbose:
        print(lengths, cents)
    return lengths, cents


def load_parameters(filename):

    with open(filename, "r") as f:
        traj = json.load(f)
    return traj


def load_3D_simus(folder_roots, n=5, S=False, skip=[], single=False, orip=False):
    from replication.ensembleSim import ensembleSim

    found = False
    if single:
        parameters = load_parameters(folder_roots + "/params.json")
        found = True
    else:
        for i in range(1, n + 1):
            try:
                parameters = load_parameters(folder_roots + "%i/params.json" % i)
                found = True
            except:
                pass
    if not found:
        print("Could not find any parameter")
        raise
    # pprint.pprint(parameters)
    if type(parameters["len_chrom"]) == str:
        lengths, _ = load_lengths_and_centro(
            "../." + parameters["len_chrom"], parameters["coarse"], verbose=False)
    else:
        lengths = parameters["len_chrom"]

    fact = 1
    if not parameters["diff_bind_when_free"]:
        fact = 2
        print("Doubling number of diff")
    E = ensembleSim(n, Nori=None,
                    Ndiff=parameters["N_diffu"] * fact,
                    lengths=lengths,
                    p_on=parameters["p_inte"],
                    p_off=parameters["p_off"],
                    only_one=parameters["diff_bind_when_free"],
                    all_same_ori=True,
                    fork_speed=parameters["fork_speed"],
                    dt_speed=parameters["dt_speed"])
    s = E.run_all(load_from_file=folder_roots, skip=skip, single=single, orip=orip)
    if S:
        return E, lengths, parameters, s

    return E, lengths, parameters
