import pandas


def load_ori_position(parameters):

    Oris = pandas.read_csv(parameters["Nori"], comment="#")

    l_ori = [[] for i in range(16)]

    # Filters ori:
    for ch, p, status in zip(Oris["chr"], Oris["start"], Oris["status"]):
        if status in parameters["ori_type"]:  # ,"Likely"]:
            l_ori[ch - 1].append(int(p / (parameters["coarse"] * 1000)))

    # Then remove duplicate (/kb ) and outside of boundaries

    tot = 0
    for i in range(len(parameters["lengths"])):
        isize = len(l_ori[i])
        l_ori[i] = list(set(l_ori[i]))
        l_ori[i].sort()
        print(isize, len(l_ori[i]))
        while not (max(l_ori[i]) < parameters["lengths"][i]):
            l_ori[i].remove(max(l_ori[i]))
        # print(max(l_ori[i]),len_chrom[i]*5

        tot += len(l_ori[i])
        assert(max(l_ori[i]) < parameters["lengths"][i])

    return l_ori
