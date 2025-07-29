import sys

def main():
    insplit_tab = sys.argv[1] # AA.ori.ccs.asm.split.fq.melTranscript.noTE.melTE.win.uniq.s.paf.combine.tab
    hap1_split_tab = sys.argv[2] # /rd/caiya/LTR/raw_ccs/AA.ori.ccs.fq.AAC.hap1.win.s.uniq.paf.asm.combine.shared.chimeric.split.tab

    out = insplit_tab.split("/")[-1] + ".asm_combine.asm_pos.tab"

    oriccs_hap1_split_dic = read_hap1_split_tab(hap1_split_tab)
    split_tab_info_dic = read_split_tab(insplit_tab)
    combine_split_tab(oriccs_hap1_split_dic, split_tab_info_dic, out)

def read_hap1_split_tab(insplit_tab):
    oriccs_dic = {}
    with open(insplit_tab, "r") as inf:
        for line in inf:
            line = line.rstrip()
            ll = line.split("\t")
            rname = ll[0]
            oriccs_dic[rname] = ll
    return oriccs_dic

# def read_hap1_split_read_bed(hap1_split_read_bed):
#     with open(hap1_split_read_bed, "r") as inf:
#         oriccs_bed_dic = {}
#         for line in inf:
#             line = line.rstrip()
#             ll = line.split("\t")
#             rname = ll[0]
#             sep_name = int(ll[3].split("_")[-1])
#             query_pos_range = ll[1] + "\t" + ll[2]
            
#             oriccs_bed_dic.setdefault(rname, {})[sep_name] = query_pos_range
#     return oriccs_bed_dic

def read_split_tab(insplit_tab):
    split_tab_info_dic = {}
    with open(insplit_tab, "r") as inf:
        for line in inf:
            line = line.rstrip()
            ll = line.split("\t")

            whole_name = "_".join(ll[0].split("_")[:-1])
            part_name = int(ll[0].split("_")[-1])

            split_tab_info_dic.setdefault(whole_name, {})[part_name] = ll

    return split_tab_info_dic


def combine_split_tab(oriccs_dic, split_tab_info_dic, out):
    with open(out, "w") as outf:
        for whole_name in oriccs_dic:
            if whole_name in split_tab_info_dic:
                oriccs_ll = oriccs_dic[whole_name]

                out_name = whole_name
                out_length_total = oriccs_ll[1]
                out_query_rnage_ls = "|".join(oriccs_ll[2].split(","))
                out_query_length_ls = "|".join(oriccs_ll[3].split(","))

                # out_query_rnage_ls = []
                # out_query_length_ls = []
                out_target_ls = []
                out_target_range_ls = []
                out_target_length_ls = []
                out_strand_ls = []
                out_type_ls = []
                out_target_name_ls = []
                out_target_type_ls = []

                for i in range(len(oriccs_ll[2].split(","))):
                    partname = i + 1
                    if partname in split_tab_info_dic[whole_name]:
                        part_ll = split_tab_info_dic[whole_name][partname]

                        out_target_ls.append(part_ll[4])
                        out_target_range_ls.append(part_ll[5])
                        out_target_length_ls.append(part_ll[6])
                        out_strand_ls.append(part_ll[7])
                        out_type_ls.append(part_ll[8])
                        out_target_name_ls.append(part_ll[9])
                        out_target_type_ls.append(part_ll[10])
                    else:
                        out_target_ls.append("unknown")
                        out_target_range_ls.append("unknown")
                        out_target_length_ls.append("unknown")
                        out_strand_ls.append("unknown")
                        out_type_ls.append("unknown")
                        out_target_name_ls.append("unknown")
                        out_target_type_ls.append("unknown")

                out_line = out_name + "\t" + out_length_total + "\t" + out_query_rnage_ls + "\t" + out_query_length_ls + "\t" + "|".join(out_target_ls) + "\t" + "|".join(out_target_range_ls) + "\t" + "|".join(out_target_length_ls) + "\t" + "|".join(out_strand_ls) + "\t" + "|".join(out_type_ls) + "\t" + "|".join(out_target_name_ls) + "\t" + "|".join(out_target_type_ls)
                outf.write(out_line + "\n")


if __name__ == "__main__":
    main()


