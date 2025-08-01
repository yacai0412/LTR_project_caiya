import sys

def main():
    ltr_type_snv = sys.argv[1] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/ecc/LTR_type/AA.deepconsensus.cx3.q30.ecc1LTR.chimeric.split.AAC.F904.shared.hap1.ds.plus_minus.subpass_5.vaf_0.8.rm_25bp.LTR3.breakpoint.BEL.snv

    out1 = ".".join(ltr_type_snv.split("/")[-1].split(".")[:-1]) + ".LTR3.breakpoint.position.out"
    breakpoint_position_dis_dic = initial_breakpoint_position_dis_dic()
    breakpoint_position_dis_dic = get_mutation_relative_pos(ltr_type_snv, breakpoint_position_dis_dic)
    write_breakpoint_position_dis_dic(breakpoint_position_dis_dic, out1)


def initial_breakpoint_position_dis_dic():
    breakpoint_position_dis_dic = {}
    for i in range(10000):
        breakpoint_position_dis_dic.setdefault("LTR_3", {})[-1 * (i + 1)] = 0
        breakpoint_position_dis_dic.setdefault("I_5", {})[i + 1] = 0
    return breakpoint_position_dis_dic


def get_mutation_relative_pos(ltr_type_snv, breakpoint_position_dis_dic):
    with open(ltr_type_snv, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            ltr_type = ll[4]
            ltr_brakpoint_pos = int(ll[5])
            breakpoint_position_dis_dic[ltr_type][ltr_brakpoint_pos] += 1
    return breakpoint_position_dis_dic

def write_breakpoint_position_dis_dic(breakpoint_position_dis_dic, out1):
    with open(out1, "w") as final:
        ltr_3_dic = breakpoint_position_dis_dic["LTR_3"]
        i_5_dic = breakpoint_position_dis_dic["I_5"]
        for i in ltr_3_dic:
            out_line = "LTR_3\t" + str(i) + "\t" + str(breakpoint_position_dis_dic["LTR_3"][i])
            final.write(out_line + "\n")

        for j in i_5_dic:
            out_line = "I_5\t" + str(j) + "\t" + str(breakpoint_position_dis_dic["I_5"][j])
            final.write(out_line + "\n")

if __name__ == "__main__":
    main()
