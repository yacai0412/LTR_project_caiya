import sys

def main():
    mutation_type_list = sys.argv[1] # AA.LR200.mutationtype
    depth = sys.argv[2] # AA.deepconsensus.cx3.dm6.F904.overlap90.q20_q30.depth

    out = mutation_type_list + ".rate"

    depth_dic = read_depth(depth)
    ds_snv_dic, ss_snv_dic = read_mutation_type_list(mutation_type_list)
    calculate_mutaetion_rate(ds_snv_dic, ss_snv_dic, depth_dic, out)


def read_depth(depth):
    depth_dic = {}
    with open(depth, "r") as inf:
        for line in inf:
            ll = line.rstrip("\n").split("\t")
            qv = ll[1]
            depth_dic.setdefault(qv, {})[ll[2]] = int(ll[3])
    return depth_dic

def read_mutation_type_list(mutation_type_list):
    ds_snv_dic = {}
    ss_snv_dic = {}
    with open(mutation_type_list, "r") as inf:
        for line in inf:
            ll = line.rstrip("\n").split("\t")
            qv = ll[1]
            if ll[0] == "ds":
                ds_snv_dic.setdefault(qv, {})[ll[2]] = int(ll[3])
            elif ll[0] == "ss":
                ss_snv_dic.setdefault(qv, {})[ll[2]] = int(ll[3])
    return ds_snv_dic, ss_snv_dic

def calculate_mutaetion_rate(ds_snv_dic, ss_snv_dic, depth_dic, out):
    type_list = ["LTR", "LINE", "DNA", "RC", "RNA", "VNTR", "Satellite", "singleton", "multicopy", "UTR", "small_8_30_intron", "other_intron", "intergenic", "Other", "Unknown"]
    
    with open(out, "w") as final:
        for i in ["q20", "q30"]:
            for j in type_list:
                if i in ds_snv_dic and i in depth_dic:
                    if j in ds_snv_dic[i] and j in depth_dic[i]:
                        snvcount = ds_snv_dic[i][j]
                        depth = depth_dic[i][j]
                        mutation_rate = snvcount / (depth / 2)
                        out_line = "ds\t" + i + "\t" + j + "\t" + str("{:.2e}".format(snvcount)) + "\t" + str("{:.2e}".format(depth)) + "\t" + str("{:.2e}".format(mutation_rate))
                        final.write(out_line + "\n")
                    else:
                        snvcount = ds_snv_dic[i].get(j, 0)
                        depth = depth_dic[i].get(j, 0)
                        mutation_rate = 0
                        out_line = "ds\t" + i + "\t" + j + "\t" + str("{:.2e}".format(snvcount)) + "\t" + str("{:.2e}".format(depth)) + "\t" + str("{:.2e}".format(mutation_rate))
                        final.write(out_line + "\n")
            final.write("\n")

        for i in ["q20", "q30"]:
            for j in type_list:
                if i in ss_snv_dic and i in depth_dic:
                    if j in ss_snv_dic[i] and j in depth_dic[i]:
                        snvcount = ss_snv_dic[i][j]
                        depth = depth_dic[i][j]
                        mutation_rate = snvcount / depth
                        out_line = "ss\t" + i + "\t" + j + "\t" + str("{:.2e}".format(snvcount)) + "\t" + str("{:.2e}".format(depth)) + "\t" + str("{:.2e}".format(mutation_rate))
                        final.write(out_line + "\n")
                    else:
                        snvcount = ss_snv_dic[i].get(j, 0)
                        depth = depth_dic[i].get(j, 0)
                        mutation_rate = 0
                        out_line = "ss\t" + i + "\t" + j + "\t" + str("{:.2e}".format(snvcount)) + "\t" + str("{:.2e}".format(depth)) + "\t" + str("{:.2e}".format(mutation_rate))
                        final.write(out_line + "\n")
            final.write("\n")


if __name__ == "__main__":
    main()
