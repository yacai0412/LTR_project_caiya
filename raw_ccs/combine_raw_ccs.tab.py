import sys

def main():
    insplit_tab = sys.argv[1] # AA.ori.ccs.asm.split.fq.melTranscript.noTE.melTE.win.uniq.s.paf.combine.tab
    hap1_split_tab = sys.argv[2] # /rd/caiya/LTR/raw_ccs/AA.ori.ccs.fq.AAC.hap1.win.s.uniq.paf.asm.combine.shared.chimeric.split.tab
    hap1_split_read_bed = sys.argv[3] # /rd/caiya/LTR/raw_ccs/AA.ori.ccs.fq.AAC.hap1.win.s.uniq.paf.asm.combine.shared.chimeric.split.bed

    out = insplit_tab.split("/")[-1] + ".asm_combine.tab"

    oriccs_hap1_split_dic = read_hap1_split_tab(hap1_split_tab)
    oriccs_hap1_split_bed_dic = read_hap1_split_read_bed(hap1_split_read_bed)
    split_tab_info_dic = read_split_tab(insplit_tab)
    combine_split_tab(oriccs_hap1_split_dic, oriccs_hap1_split_bed_dic, split_tab_info_dic, out)


def read_hap1_split_tab(insplit_tab):
    oriccs_dic = {}
    with open(insplit_tab, "r") as inf:
        for line in inf:
            line = line.rstrip()
            ll = line.split("\t")
            rname = ll[0]
            oriccs_dic[rname] = ll
    return oriccs_dic

def read_hap1_split_read_bed(hap1_split_read_bed):
    with open(hap1_split_read_bed, "r") as inf:
        oriccs_bed_dic = {}
        for line in inf:
            line = line.rstrip()
            ll = line.split("\t")
            rname = ll[0]
            sep_name = int(ll[3].split("_")[-1])
            query_pos_range = ll[1] + "\t" + ll[2]
            
            oriccs_bed_dic.setdefault(rname, {})[sep_name] = query_pos_range
    return oriccs_bed_dic


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


def combine_split_tab(oriccs_dic, oriccs_bed_dic, split_tab_info_dic, out):
    with open(out, "w") as outf:
        for whole_name in split_tab_info_dic:
            if whole_name in oriccs_dic and whole_name in oriccs_bed_dic:
                oriccs_ll = oriccs_dic[whole_name]
        
                out_name = whole_name
                out_length_total = oriccs_ll[1]
                out_query_rnage_ls = []
                out_query_length_ls = []
                out_target_ls = []
                out_target_range_ls = []
                out_target_length_ls = []
                out_strand_ls = []
                out_type_ls = []
                out_target_name_ls = []
                out_target_type_ls = []

                for part_name in sorted(split_tab_info_dic[whole_name]):
                    ll = split_tab_info_dic[whole_name][part_name]
                    # query_startpos = int(oriccs_ll[2].split(",")[part_name - 1].split("-")[0])
                    query_startpos = int(oriccs_bed_dic[whole_name][part_name].split("\t")[0])

                    query_rnage_ls = []
                    for i in ll[2].split(","):
                        query_range_start_i = int(i.split("-")[0]) + query_startpos
                        query_range_end_i = int(i.split("-")[1]) + query_startpos
                        query_rnage_i = str(query_range_start_i) + "-" + str(query_range_end_i)
                        query_rnage_ls.append(query_rnage_i)

                    query_rnage = ",".join(query_rnage_ls)
                    out_query_rnage_ls.append(query_rnage)

                    out_query_length_ls.append(ll[3])
                    out_target_ls.append(ll[4])
                    out_target_range_ls.append(ll[5])
                    out_target_length_ls.append(ll[6])
                    out_strand_ls.append(ll[7])
                    out_type_ls.append(ll[8])
                    out_target_name_ls.append(ll[9])
                    out_target_type_ls.append(ll[10])

                out_line = out_name + "\t" + out_length_total + "\t" + "|".join(out_query_rnage_ls) + "\t" + "|".join(out_query_length_ls) + "\t" + "|".join(out_target_ls) + "\t" + "|".join(out_target_range_ls) + "\t" + "|".join(out_target_length_ls) + "\t" + "|".join(out_strand_ls) + "\t" + "|".join(out_type_ls) + "\t" + "|".join(out_target_name_ls) + "\t" + "|".join(out_target_type_ls)
                outf.write(out_line + "\n")

if __name__ == "__main__":
    main()
