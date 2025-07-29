import sys
import re

def main():
    in_tab1 = sys.argv[1] # AAC.ccs.deepconsensus.asm.split.fq.melTE.mm2sr.s.paf.combine.tab1

    out_vcf = in_tab1 + ".TEI.vcf"

    te_ins_dic = read_intab1_get_tEI(in_tab1)
    export_tei_ins_vcf(te_ins_dic, out_vcf)


def normalize_target_name(name):
    if 'HMSBEAGLE' in name or 'DMLTR5' in name:
        return 'HMSBEAGLE'
    else:
        return name.split('_LTR')[0].split("-LTR")[0].split("-I")[0].split("_I")[0]  # Extract the base name before '_'


def read_intab1_get_tEI(in_tab1):
    te_ins_dic = {}

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
                            
                            te_ins_dic.setdefault(ins_pos, []).append((rname, ltr_names[0]))


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
                            
                            te_ins_dic.setdefault(ins_pos, []).append((rname, ltr_names[0]))


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

                            te_ins_dic.setdefault(ins_pos, []).append((rname, ltr_names[0]))

    return te_ins_dic


def export_tei_ins_vcf(te_ins_dic, out_vcf):
    # Prepare list of (chrom, pos, ltr_family, list_of_reads)
    tei_list = []
    for ins_pos in te_ins_dic:
        chrom, pos = ins_pos.split("\t")
        pos = int(pos)
        out_put_dic = {}
        for rname, ltr_family in te_ins_dic[ins_pos]:
            out_put_dic.setdefault(ltr_family, []).append(rname)
        for ltr_family in out_put_dic:
            tei_list.append( (chrom, pos, ltr_family, out_put_dic[ltr_family]) )

    # Sort by chrom and position
    tei_list.sort(key=lambda x: (x[0], x[1]))

    # Merge close insertions
    merged_list = []
    if tei_list:
        last_chrom, last_pos, last_ltr_family, last_reads = tei_list[0]
        for i in range(1, len(tei_list)):
            chrom, pos, ltr_family, reads = tei_list[i]
            if chrom == last_chrom and ltr_family == last_ltr_family and abs(pos - last_pos) <= 10:
                last_reads.extend(reads)
                last_pos = min(last_pos, pos)  # always keep smallest pos
            else:
                merged_list.append( (last_chrom, last_pos, last_ltr_family, sorted(set(last_reads))) )
                last_chrom, last_pos, last_ltr_family, last_reads = chrom, pos, ltr_family, reads
        # Add last
        merged_list.append( (last_chrom, last_pos, last_ltr_family, sorted(set(last_reads))) )

    # Write output
    with open(out_vcf, "w") as outf:
        for chrom, pos, ltr_family, reads in merged_list:
            outline = f"{chrom}\t{pos}\t{ltr_family}\t{','.join(reads)}\n"
            outf.write(outline)


if __name__ == "__main__":
    main()
