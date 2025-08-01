import sys

def main():
    septab = sys.argv[1] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/ecc/AA.deepconsensus.cx3.fwd_rev.trans1_te1_ref.win.s.paf.ecc1LTR.rm_no_split.tab
    snv = sys.argv[2] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/ecc/AA.deepconsensus.cx3.q30.ecc1LTR.chimeric.split.AAC.F904.shared.hap1.ds.plus_minus.subpass_5.vaf_0.8.rm_25bp.out
    ds_ss_type = sys.argv[3] # ds

    out = ".".join(snv.split("/")[-1].split(".")[:-1]) + ".LTR3.breakpoint.snv"
    out1 = ".".join(snv.split("/")[-1].split(".")[:-1]) + ".LTR3.breakpoint.position.out"


    breakpoint_position_dis_dic = initial_breakpoint_position_dis_dic()

    sep_1LTR_plus_part_length_dic, sep_1LTR_minus_part_length_dic = read_tab_get_plus_sep_LTR(septab)
    if ds_ss_type == "ds":
        breakpoint_position_dis_dic = read_ds_snv(snv, sep_1LTR_plus_part_length_dic, breakpoint_position_dis_dic, out)
    elif ds_ss_type == "ss":
        breakpoint_position_dis_dic = read_ss_snv(snv, sep_1LTR_plus_part_length_dic, sep_1LTR_minus_part_length_dic, breakpoint_position_dis_dic, out)
    write_breakpoint_position_dis_dic(breakpoint_position_dis_dic, out1)


def read_tab_get_plus_sep_LTR(septab):
    sep_1LTR_plus_part_length_dic = {}
    sep_1LTR_minus_part_length_dic = {}
    with open(septab, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            rname = ll[0]
            part_length = [int(x) for x in ll[3].split(",")]
            te_type = ll[4]

            if ll[7] == "+,+,+":
                sep_1LTR_plus_part_length_dic[rname] = part_length, te_type
            else:
                sep_1LTR_minus_part_length_dic[rname] = part_length, te_type
    return sep_1LTR_plus_part_length_dic, sep_1LTR_minus_part_length_dic

def initial_breakpoint_position_dis_dic():
    breakpoint_position_dis_dic = {}
    for i in range(10000):
        breakpoint_position_dis_dic.setdefault("LTR_3", {})[-1 * (i + 1)] = 0
        breakpoint_position_dis_dic.setdefault("I_5", {})[i + 1] = 0
    return breakpoint_position_dis_dic

def read_ds_snv(snv, sep_1LTR_plus_part_length_dic, breakpoint_position_dis_dic, out):
    with open(snv, "r") as inf, open(out, "w") as final:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            out_line_part0 = "\t".join(ll[0:4])
            plus_reads_ls = ll[6].split(",")
            minus_reads_ls = ll[7].split(",")

            for i in range(len(plus_reads_ls)):
                plus_rname0 = plus_reads_ls[i]
                minus_rname0 = minus_reads_ls[i]
                
                plus_rname1 = "_".join(plus_rname0.split(":")[0].split("_")[:-1])
                minus_rname1 = "_".join(minus_rname0.split(":")[0].split("_")[:-1])

                if plus_rname1 in sep_1LTR_plus_part_length_dic:
                    ltr_plus_rname = plus_rname0
                    part_readlength_ls = sep_1LTR_plus_part_length_dic[plus_rname1][0]
                    te_type = sep_1LTR_plus_part_length_dic[plus_rname1][1]
                elif minus_rname1 in sep_1LTR_plus_part_length_dic:
                    ltr_plus_rname = minus_rname0
                    part_readlength_ls = sep_1LTR_plus_part_length_dic[minus_rname1][0]
                    te_type = sep_1LTR_plus_part_length_dic[minus_rname1][1]
                else:
                    print("error in " + line + "\t" + plus_rname1)
                    continue

                sep_part = int(ltr_plus_rname.split(":")[0].split("_")[-1])
                sep_pos = int(ltr_plus_rname.split(":")[1])

                if sep_part == 2:
                    sep_ltr_type = "LTR_3"
                    part_length = part_readlength_ls[1]
                    breakpoint_distance = sep_pos - part_length
                    breakpoint_position_dis_dic[sep_ltr_type][breakpoint_distance] += 1
                    out_line = out_line_part0 + "\t" + sep_ltr_type + "\t" + str(breakpoint_distance) + "\t" + te_type + "\t" + plus_rname0 + "," + minus_rname0
                    final.write(out_line + "\n")
                elif sep_part == 3:
                    sep_ltr_type = "I_5"
                    breakpoint_distance = sep_pos
                    breakpoint_position_dis_dic[sep_ltr_type][breakpoint_distance] += 1
                    out_line = out_line_part0 + "\t" + sep_ltr_type + "\t" + str(breakpoint_distance) + "\t" + te_type + "\t" + plus_rname0 + "," + minus_rname0
                    final.write(out_line + "\n")

    return breakpoint_position_dis_dic


def read_ss_snv(snv, sep_1LTR_plus_part_length_dic, sep_1LTR_minus_part_length_dic, breakpoint_position_dis_dic, out):
    with open(snv, "r") as inf, open(out, "w") as final:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            out_line_part0 = "\t".join(ll[0:4])
            support_reads_ls = ll[6].split(",")

            for i in range(len(support_reads_ls)):
                support_rname0 = support_reads_ls[i]
                support_rname1 = "_".join(support_rname0.split(":")[0].split("_")[:-1])

                sep_part = int(support_rname0.split(":")[0].split("_")[-1])
                sep_pos = int(support_rname0.split(":")[1])

                if support_rname1 in sep_1LTR_plus_part_length_dic:
                    part_readlength_ls = sep_1LTR_plus_part_length_dic[support_rname1][0]
                    te_type = sep_1LTR_plus_part_length_dic[support_rname1][1]
                    if sep_part == 2:
                        sep_ltr_type = "LTR_3"
                        part_length = part_readlength_ls[1]
                        breakpoint_distance = sep_pos - part_length
                        breakpoint_position_dis_dic[sep_ltr_type][breakpoint_distance] += 1
                        out_line = out_line_part0 + "\t" + sep_ltr_type + "\t" + str(breakpoint_distance) + "\t" + te_type + "\t" + support_rname0
                        final.write(out_line + "\n")
                    elif sep_part == 3:
                        sep_ltr_type = "I_5"
                        breakpoint_distance = sep_pos
                        breakpoint_position_dis_dic[sep_ltr_type][breakpoint_distance] += 1
                        out_line = out_line_part0 + "\t" + sep_ltr_type + "\t" + str(breakpoint_distance) + "\t" + te_type + "\t" + support_rname0
                        final.write(out_line + "\n")

                elif support_rname1 in sep_1LTR_minus_part_length_dic:
                    part_readlength_ls = sep_1LTR_minus_part_length_dic[support_rname1][0]
                    te_type = sep_1LTR_minus_part_length_dic[support_rname1][1]
                    if sep_part == 2:
                        sep_ltr_type = "LTR_3"
                        breakpoint_distance = sep_pos * (-1)
                        breakpoint_position_dis_dic[sep_ltr_type][breakpoint_distance] += 1
                        out_line = out_line_part0 + "\t" + sep_ltr_type + "\t" + str(breakpoint_distance) + "\t" + te_type + "\t" + support_rname0
                        final.write(out_line + "\n")
                    elif sep_part == 1:
                        sep_ltr_type = "I_5"
                        part_length = part_readlength_ls[0]
                        breakpoint_distance = part_length - sep_pos
                        breakpoint_position_dis_dic[sep_ltr_type][breakpoint_distance] += 1
                        out_line = out_line_part0 + "\t" + sep_ltr_type + "\t" + str(breakpoint_distance) + "\t" + te_type + "\t" + support_rname0
                        final.write(out_line + "\n")
                else:
                    print("error in " + support_rname1 + " can not find in tab")


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







