import pandas


def load_ori_position(File, ori_type, lengths, coarse):

    Oris = pandas.read_csv(File, comment="#")

    l_ori = [[] for i in range(len(lengths))]

    # Filters ori:
    for ch, p, status in zip(Oris["chr"], Oris["start"], Oris["status"]):
        if status in ori_type:  # ,"Likely"]:
            l_ori[ch - 1].append(int(p / coarse))

    # Then remove duplicate (/kb ) and outside of boundaries

    tot = 0
    for i in range(len(lengths)):
        isize = len(l_ori[i])
        l_ori[i] = list(set(l_ori[i]))
        l_ori[i].sort()
        print(isize, len(l_ori[i]))
        while not (max(l_ori[i]) < lengths[i]):
            l_ori[i].remove(max(l_ori[i]))
        # print(max(l_ori[i]),len_chrom[i]*5

        tot += len(l_ori[i])
        assert(max(l_ori[i]) < lengths[i])

    return l_ori


def get_lengths_and_centro(File, coarse=1):
    gff = pandas.read_csv(File, sep="\t", comment="#", header=None, names=[
                          "chr", "SGD", "info", "start", "end", "m1", "m2", "m3", "comment"])

    lengths = [int(end) // int(coarse) for inf, end, chro in zip(gff["info"], gff["end"], gff["chr"])
               if inf == "chromosome" and "mt" not in chro]
    # for i in range(remove):
    # print(lengths,len(lengths))
    cents = [int(end) // int(coarse)
             for inf, end in zip(gff["info"], gff["end"]) if inf == "centromere"]
    # print(cents,len(cents))
    return lengths, cents
