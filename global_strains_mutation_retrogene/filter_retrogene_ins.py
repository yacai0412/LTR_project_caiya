import sys
import re

def main():
    in_gene_tei_tab1 = sys.argv[1] # G0-F100.ont.asm.split.fq.melTE.mm2sr.s.paf.combine.tab1.TEI.gene.tab1
    sep_LTR_fasta = "/public5/home/sch5655/global_strains/melTranscript.noTE.melTE.fa"

    out = in_gene_tei_tab1.replace("tab1", "tab2")

    te_ref_length_dic = read_sep_LTR_fasta_get_legnth(sep_LTR_fasta)
    read_in_gene_tei_tab1(in_gene_tei_tab1, te_ref_length_dic, out)


def read_sep_LTR_fasta_get_legnth(sep_LTR_fasta):
    te_ref_length_dic = {}
    with open(sep_LTR_fasta, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line.strip(">").split("|")[0]
                te_ref_length_dic[name] = 0
            else:
                te_ref_length_dic[name] += len(line.strip())
    return te_ref_length_dic

def normalize_target_name(name):
    if 'HMSBEAGLE' in name or 'DMLTR5' in name:
        return 'HMSBEAGLE'
    else:
        return name.split('_LTR')[0].split("-LTR")[0].split("-I")[0].split("_I")[0]  # Extract the base name before '_'


def if_LTR_complete(LTR_index, segments, target_range_ls, te_ref_length_dic):
    cutoff = 15
    ltr_name = segments[LTR_index]
    ltr_target_range = target_range_ls[LTR_index].split("-")

    if ltr_name in te_ref_length_dic:
        whole_length = te_ref_length_dic[ltr_name]

        # print(ltr_name)
        # print(ltr_target_range)
        # print(whole_length)

        if 0 <= int(ltr_target_range[0]) <= cutoff and (whole_length - cutoff) <= int(ltr_target_range[1]) <= whole_length:
            return True
        else:
            return False
    else:
        return False


def read_in_gene_tei_tab1(in_gene_tei_tab1, te_ref_length_dic, out):
    with open(in_gene_tei_tab1, "r") as inf, open(out, "w") as outf:
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
            geneseg = [i for i, seg in enumerate(segments) if 'FBgn' in seg]
            
            if (len(chroms) + len(ltrs) + len(iseg) + len(geneseg)) == len(segments):
                if len(chroms) > 0 and len(ltrs) > 0 and len(iseg) > 0:
                    ltr_names = [normalize_target_name(segments[i]) for i in ltrs]
                    i_names   = [normalize_target_name(segments[i]) for i in iseg]
                    same_ltr_family = all(name == ltr_names[0] for name in ltr_names)
                    same_i_family   = all(name == i_names[0] for name in i_names)
                    same_family     = same_ltr_family and same_i_family and (ltr_names[0] == i_names[0])

                    same_ltr_i_strand = len(set([strand_ls[i] for i in ltrs + iseg])) == 1
                    same_chr_strand   = len(set([strand_ls[i] for i in chroms])) == 1


                    # if (len(chroms) >= 2 and len(ltrs) == 2 and len(iseg) >=1):
                    if len(ltrs) == 2 and len(iseg) >= 1:
                        # Ensure order: chr - LTR - I - LTR - chr
                        first_chr_group = [i for i in chroms if i < ltrs[0]]
                        last_chr_group = [i for i in chroms if i > ltrs[-1]]

                        left_LTR = min(ltrs + iseg)
                        rigth_LTR = max(ltrs + iseg)


                        if len(first_chr_group) >= 1 and len(last_chr_group) >= 1:
                            chr_names = [segments[i] for i in first_chr_group + last_chr_group]
                            same_chr = all(name == chr_names[0] for name in chr_names)

                            ordered = first_chr_group[-1] < ltrs[0] < iseg[0] < ltrs[-1] < last_chr_group[0]
                            proximity = (first_chr_group[-1] + 1 == ltrs[0]) and (ltrs[-1] + 1 == last_chr_group[0])
                            surrounded = all(ltrs[0] < i < ltrs[-1] for i in iseg)
                            gene_inside = all(left_LTR < g < rigth_LTR for g in geneseg)
                            complete_LTR = if_LTR_complete(ltrs[0], segments, target_range_ls, te_ref_length_dic) and if_LTR_complete(ltrs[-1], segments, target_range_ls, te_ref_length_dic)

                            if ordered and proximity and surrounded and same_chr and same_family and same_ltr_i_strand and same_chr_strand and gene_inside and complete_LTR:
                                outf.write(line + "\n")


                    elif (len(chroms) >=1 and len(ltrs) >=1 and len(iseg) >= 1):
                        first_chr_group = [i for i in chroms if i < ltrs[0]]
                        last_chr_group = [i for i in chroms if i > ltrs[-1]]

                        left_LTR = min(ltrs + iseg)
                        rigth_LTR = max(ltrs + iseg)

                        chr_names = [segments[i] for i in first_chr_group + last_chr_group]
                        same_chr = all(name == chr_names[0] for name in chr_names)

                        if first_chr_group and first_chr_group[-1] + 1 == ltrs[0] and last_chr_group == []:
                            # chr-LTR-I
                            ordered = first_chr_group[-1] < ltrs[0] < iseg[0]
                            proximity = first_chr_group[-1] + 1 == ltrs[0]
                            surrounded = all(ltrs[0] < i for i in iseg)
                            gene_inside = all(left_LTR < g for g in geneseg)
                            complete_LTR = if_LTR_complete(ltrs[0], segments, target_range_ls, te_ref_length_dic)

                            if ordered and proximity and surrounded and same_chr and same_family and same_ltr_i_strand and same_chr_strand and gene_inside and complete_LTR:
                                outf.write(line + "\n")


                        elif last_chr_group and ltrs[-1] + 1 == last_chr_group[0] and first_chr_group == []:
                            # I-LTR-chr
                            ordered = iseg[0] < ltrs[-1] < last_chr_group[0]
                            proximity = ltrs[-1] + 1 == last_chr_group[0]
                            surrounded = all(i < ltrs[-1] for i in iseg)
                            gene_inside = all(g < rigth_LTR for g in geneseg)
                            complete_LTR = if_LTR_complete(ltrs[-1], segments, target_range_ls, te_ref_length_dic)

                            if ordered and proximity and surrounded and same_chr and same_family and same_ltr_i_strand and same_chr_strand and gene_inside and complete_LTR:
                                outf.write(line + "\n")

if __name__ == "__main__":
    main()

