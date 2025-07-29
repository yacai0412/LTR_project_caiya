import sys
import re


def main():
    combine_all_tab1 = sys.argv[1] # /public5/home/sch5655/global_strains/intron_loss_retrogene/G0.ont.asm.split.fq.melTE.mm2sr.s.paf.combine.tab1
    split_gene_tab = sys.argv[2] # /public5/home/sch5655/global_strains/intron_loss_retrogene/G0.asm.split.TEI.fq.melTranscript.noTE.melTE.mm2ont.s.paf.combine.gene.tab

    out = ".".join(combine_all_tab1.split("/")[-1].split(".")[:-1]) + ".TEI.gene.tab1"

    te_ins_parts_dic = read_intab1_get_tEI(combine_all_tab1)
    split_reads_parts_dic = read_split_tab(split_gene_tab)
    combine_tab1_split_gene_dic(te_ins_parts_dic, split_reads_parts_dic, out)


def normalize_target_name(name):
    if 'HMSBEAGLE' in name or 'DMLTR5' in name:
        return 'HMSBEAGLE'
    else:
        return name.split('_LTR')[0].split("-LTR")[0].split("-I")[0].split("_I")[0]  # Extract the base name before '_'

def read_intab1_get_tEI(in_tab1):
    # te_ins_dic = {}
    te_ins_parts_dic = {}

    with open(in_tab1, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            
            rname = ll[0]
            target_range_ls = re.split(r'[|,]', ll[5])
            strand_ls = re.split(r'[|,]', ll[7])

            # Split by both | and ,
            segments = re.split(r'[|,]', ll[4])

            # Remove any extra spaces
            segments = [seg.strip() for seg in segments]

            chroms = [i for i, seg in enumerate(segments) if seg.startswith('chr')]
            ltrs = [i for i, seg in enumerate(segments) if 'LTR' in seg]
            iseg = [i for i, seg in enumerate(segments) if '-I' in seg or '_I' in seg]


            if len(chroms) > 0 and len(ltrs) > 0 and len(iseg) > 0:
                ltr_names = [normalize_target_name(segments[i]) for i in ltrs]
                i_names   = [normalize_target_name(segments[i]) for i in iseg]
                same_ltr_family = all(name == ltr_names[0] for name in ltr_names)
                same_i_family   = all(name == i_names[0] for name in i_names)
                same_family     = same_ltr_family and same_i_family and (ltr_names[0] == i_names[0])

                same_ltr_i_strand = len(set([strand_ls[i] for i in ltrs + iseg])) == 1
                same_chr_strand   = len(set([strand_ls[i] for i in chroms])) == 1


                # if (len(chroms) >= 2 and len(ltrs) == 2 and len(iseg) >=1):
                if len(ltrs) == 2 and len(iseg) >=1:
                    # Ensure order: chr - LTR - I - LTR - chr
                    first_chr_group = [i for i in chroms if i < ltrs[0]]
                    last_chr_group = [i for i in chroms if i > ltrs[-1]]


                    if len(first_chr_group) >= 1 and len(last_chr_group) >= 1:
                        chr_names = [segments[i] for i in first_chr_group + last_chr_group]
                        same_chr = all(name == chr_names[0] for name in chr_names)

                        ordered = first_chr_group[-1] < ltrs[0] < iseg[0] < ltrs[-1] < last_chr_group[0]
                        proximity = (first_chr_group[-1] + 1 == ltrs[0]) and (ltrs[-1] + 1 == last_chr_group[0])
                        surrounded = all(ltrs[0] < i < ltrs[-1] for i in iseg)

                        if ordered and proximity and surrounded and same_chr and same_family and same_ltr_i_strand and same_chr_strand:
                            if strand_ls[first_chr_group[-1]] == "+":
                                ins_pos = segments[chroms[0]] + "\t" + target_range_ls[first_chr_group[-1]].split("-")[1]
                            else:
                                ins_pos = segments[chroms[0]] + "\t" + target_range_ls[last_chr_group[0]].split("-")[1]
                            
                            # te_ins_dic.setdefault(ins_pos, []).append((rname, ltr_names[0]))
                            te_ins_parts_dic[rname] = ll


                elif (len(chroms) >=1 and len(ltrs) >=1 and len(iseg) >=1):
                    first_chr_group = [i for i in chroms if i < ltrs[0]]
                    last_chr_group = [i for i in chroms if i > ltrs[-1]]

                    chr_names = [segments[i] for i in first_chr_group + last_chr_group]
                    same_chr = all(name == chr_names[0] for name in chr_names)

                    if first_chr_group and first_chr_group[-1] + 1 == ltrs[0] and last_chr_group == []:
                        # chr-LTR-I
                        ordered = first_chr_group[-1] < ltrs[0] < iseg[0]
                        proximity = first_chr_group[-1] + 1 == ltrs[0]
                        surrounded = all(ltrs[0] < i for i in iseg)

                        if ordered and proximity and surrounded and same_chr and same_family and same_ltr_i_strand and same_chr_strand:
                            if strand_ls[first_chr_group[-1]] == "+":
                                ins_pos = segments[chroms[0]] + "\t" + target_range_ls[first_chr_group[-1]].split("-")[1]
                            else:
                                ins_pos = segments[chroms[0]] + "\t" + target_range_ls[first_chr_group[-1]].split("-")[0]
                            
                            # te_ins_dic.setdefault(ins_pos, []).append((rname, ltr_names[0]))
                            te_ins_parts_dic[rname] = ll


                    elif last_chr_group and ltrs[-1] + 1 == last_chr_group[0] and first_chr_group == []:
                        # I-LTR-chr
                        ordered = iseg[0] < ltrs[-1] < last_chr_group[0]
                        proximity = ltrs[-1] + 1 == last_chr_group[0]
                        surrounded = all(i < ltrs[-1] for i in iseg)

                        if ordered and proximity and surrounded and same_chr and same_family and same_ltr_i_strand and same_chr_strand:
                            if strand_ls[last_chr_group[0]] == "+":
                                ins_pos = segments[chroms[0]] + "\t" + target_range_ls[last_chr_group[0]].split("-")[0]
                            else:
                                ins_pos = segments[chroms[0]] + "\t" + target_range_ls[last_chr_group[0]].split("-")[1]

                            # te_ins_dic.setdefault(ins_pos, []).append((rname, ltr_names[0]))
                            te_ins_parts_dic[rname] = ll

    return te_ins_parts_dic


def read_split_tab(split_gene_tab):
    split_reads_parts_dic = {}

    with open(split_gene_tab, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            rname = ll[0]

            whole_name = rname.split("_")[0]  # Extract the base name before '_'
            part_name = int(rname.split("_")[-1])  # Extract the base name before '_'

            split_reads_parts_dic.setdefault(whole_name, {})[part_name] = ll
    return split_reads_parts_dic


def combine_tab1_split_gene_dic(te_ins_parts_dic, split_reads_parts_dic, out):
    with open(out, "w") as outf:
        for whole_name in te_ins_parts_dic:
            if whole_name in split_reads_parts_dic:
                whole_ranme = whole_name
                whole_length = te_ins_parts_dic[whole_name][1]
                whole_query_range_ls = te_ins_parts_dic[whole_name][2].split(",")
                whole_query_length_ls = te_ins_parts_dic[whole_name][3].split(",")
                whole_target_names_ls = te_ins_parts_dic[whole_name][4].split("|")
                whole_target_range_ls = te_ins_parts_dic[whole_name][5].split("|")
                whole_target_length_ls = te_ins_parts_dic[whole_name][6].split("|")
                whole_target_strand_ls = te_ins_parts_dic[whole_name][7].split("|")
                
                out_query_range_ls = []
                out_query_length_ls = []
                out_target_names_ls = []
                out_target_range_ls = []
                out_target_length_ls = []
                out_target_strand_ls = []

                for i in range(len(whole_query_range_ls)):
                    if (i + 1) in split_reads_parts_dic[whole_name]:
                        split_read_parts_ll = split_reads_parts_dic[whole_name][i + 1]
                        out_query_range_ls.append(split_read_parts_ll[2])
                        out_query_length_ls.append(split_read_parts_ll[3])
                        out_target_names_ls.append(split_read_parts_ll[4])
                        out_target_range_ls.append(split_read_parts_ll[5])
                        out_target_length_ls.append(split_read_parts_ll[6])
                        out_target_strand_ls.append(split_read_parts_ll[7])
                    else:
                        out_query_range_ls.append(whole_query_range_ls[i])
                        out_query_length_ls.append(whole_query_length_ls[i])
                        out_target_names_ls.append(whole_target_names_ls[i])
                        out_target_range_ls.append(whole_target_range_ls[i])
                        out_target_length_ls.append(whole_target_length_ls[i])
                        out_target_strand_ls.append(whole_target_strand_ls[i])

                out_line = f"{whole_ranme}\t{whole_length}\t{'|'.join(out_query_range_ls)}\t{'|'.join(out_query_length_ls)}\t{'|'.join(out_target_names_ls)}\t{'|'.join(out_target_range_ls)}\t{'|'.join(out_target_length_ls)}\t{'|'.join(out_target_strand_ls)}\n"
                outf.write(out_line)

if __name__ == "__main__":
    main()
