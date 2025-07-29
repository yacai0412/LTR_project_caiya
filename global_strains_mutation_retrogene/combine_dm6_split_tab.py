import sys


def main():
    whole_read_dm6_paf_split_tab = sys.argv[1] # AAC.ccs.deepconsensus.fq.dm6.mm2hifi.s.paf.asm.combine.tab
    split_read_melTE_paf_tab = sys.argv[2] # AAC.ccs.deepconsensus.asm.split.fq.melTranscript.canonical_TE.mm2sr.s.paf.combine.tab

    whole_reads_dm6_parts_dic, whole_reads_dm6_parts_num_dic = read_whole_read_dm6_paf_split_tab(whole_read_dm6_paf_split_tab)
    split_reads_parts_dic = read_split_read_melTE_paf_tab(split_read_melTE_paf_tab)

    out = split_read_melTE_paf_tab + "1"
    combine_whole_split_reads_parts_tab(whole_reads_dm6_parts_dic, whole_reads_dm6_parts_num_dic, split_reads_parts_dic, out)


def read_whole_read_dm6_paf_split_tab(whole_read_dm6_paf_split_tab):
    whole_reads_dm6_parts_dic = {}
    whole_reads_dm6_parts_num_dic = {}
    with open(whole_read_dm6_paf_split_tab, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")

            rname = ll[0] # m64311e_220111_083835/100007993/ccs
            readlength = int(ll[1]) # 10224
            query_range_ls = ll[2].split(",")

            for i in range(len(query_range_ls)):
                subread_name = f"{rname}_{i+1}"
                whole_reads_dm6_parts_dic.setdefault(rname, {})[subread_name] = ll
            whole_reads_dm6_parts_num_dic[rname] = len(query_range_ls)
    return whole_reads_dm6_parts_dic, whole_reads_dm6_parts_num_dic

def read_split_read_melTE_paf_tab(split_read_melTE_paf_tab):
    split_reads_parts_dic = {}
    with open(split_read_melTE_paf_tab, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")

            split_name = ll[0] # m64311e_220111_083835/100007993/ccs_1
            whole_name = "_".join(split_name.split("_")[:-1])

            split_reads_parts_dic.setdefault(whole_name, {})[split_name] = ll
    return split_reads_parts_dic


def combine_whole_split_reads_parts_tab(whole_reads_dm6_parts_dic, whole_reads_dm6_parts_num_dic, split_reads_parts_dic, out):
    with open(out, "w") as outf:
        for whole_name in whole_reads_dm6_parts_num_dic:
            whole_reads_dm6_parts_num = whole_reads_dm6_parts_num_dic[whole_name]

            whole_read_parts_ll = whole_reads_dm6_parts_dic[whole_name][whole_name + "_1"]
            out_rname = whole_name
            out_readlength = whole_read_parts_ll[1]
            out_query_range_ls_str = whole_read_parts_ll[2]
            out_query_len_ls_str = whole_read_parts_ll[3]
            out_target_name_ls = []
            out_target_range_ls = []
            out_target_length_ls = []
            out_strand_ls = []

            for i in range(whole_reads_dm6_parts_num):
                i_1 = i + 1
                split_read_name = whole_name + "_" + str(i_1)

                if whole_name in whole_reads_dm6_parts_dic and split_read_name in whole_reads_dm6_parts_dic[whole_name]:
                    whole_read_parts_ll = whole_reads_dm6_parts_dic[whole_name][split_read_name]

                    if whole_name in split_reads_parts_dic and split_read_name in split_reads_parts_dic[whole_name]:
                        split_read_parts_ll = split_reads_parts_dic[whole_name][split_read_name]
                        out_target_name_ls.append(split_read_parts_ll[4])
                        out_target_range_ls.append(split_read_parts_ll[5])
                        out_target_length_ls.append(split_read_parts_ll[6])
                        out_strand_ls.append(split_read_parts_ll[7])
                    else:
                        out_target_name_ls.append(whole_read_parts_ll[4].split(",")[i])
                        out_target_range_ls.append(whole_read_parts_ll[5].split(",")[i])
                        out_target_length_ls.append(whole_read_parts_ll[6].split(",")[i])
                        out_strand_ls.append(whole_read_parts_ll[7].split(",")[i])
                else:
                    print("error in " + whole_name)
            
            out_line = "\t".join([
                out_rname,
                str(out_readlength),
                out_query_range_ls_str,
                out_query_len_ls_str,
                "|".join(out_target_name_ls),
                "|".join(out_target_range_ls),
                "|".join(out_target_length_ls),
                "|".join(out_strand_ls),
            ]) + "\n"
            outf.write(out_line)


if __name__ == "__main__":
    main()
