import sys

def main():
    fixed_bam_LTR_bed = sys.argv[1]  # AAC_AUC.ccs.deepconsensus.asm.split.fq.melTE.mm2sr.s.paf.combine.tab1.TEI.AAC.specific.vcf.fixed.asm.bed
    hap1_rmsk = sys.argv[2]
    hap2_rmsk = sys.argv[3]

    out = ".".join(fixed_bam_LTR_bed.split("/")[-1].split(".")[:-1]) + ".overlap_rmsk.bed"

    hap12_rmsk_LTR_dic = {}
    hap12_rmsk_LTR_dic = read_rmsk_get_LTR(hap1_rmsk, hap12_rmsk_LTR_dic)
    hap12_rmsk_LTR_dic = read_rmsk_get_LTR(hap2_rmsk, hap12_rmsk_LTR_dic)
    read_fixed_bam_LTR_bed(fixed_bam_LTR_bed, hap12_rmsk_LTR_dic, out)

def read_rmsk_get_LTR(hap_rmsk, rmsk_LTR_dic):
    c = 0
    with open(hap_rmsk, "r") as inf:
        for line in inf:
            line = line.strip()
            ll = line.split()

            chrom = ll[4]
            start = int(ll[5])
            end = int(ll[6])
            strand = ll[8]
            te_name = ll[9]
            te_family = ll[10]

            if te_family.split("/")[0] == "LTR" and "LTR" in te_name:
                c += 1
                rmsk_LTR_dic.setdefault(chrom, {})[c] = (start, end, te_name, te_family, strand)
    return rmsk_LTR_dic

def is_overlap(start1, end1, start2, end2):
    return not (end1 <= start2 or end2 <= start1)

def if_same_TEname(name1, name2):
    if name1 == name2:
        return True
    elif name1 in name2 or name2 in name1:
        return True
    elif "Gypsy" in name1 and "Gypsy" in name2:
        return True
    elif ("HMS" in name1 or "DMLTR5" in name1) and ("HMS" in name2 or "DMLTR5" in name2):
        return True
    else:
        return False


def read_fixed_bam_LTR_bed(fixed_bam_LTR_bed, hap12_rmsk_LTR_dic, out):
    with open(fixed_bam_LTR_bed, "r") as inf, open(out, "w") as outf:
        for line in inf:
            line = line.strip()
            ll = line.split("\t")
            chrom = ll[0]
            start = int(ll[1])
            end = int(ll[2])
            te_name = ll[3].split(":")[2]

            fixed_LTR_bed = {}
            fixed_LTR_pos_ls = []
            if chrom in hap12_rmsk_LTR_dic:
                for c, (r_start, r_end, r_te_name, r_te_family, r_strand) in hap12_rmsk_LTR_dic[chrom].items():
                    if is_overlap(start, end, r_start, r_end) and if_same_TEname(te_name, r_te_name):
                        fixed_LTR_bed[c] = r_start, r_end, r_te_name, r_te_family, r_strand
                        fixed_LTR_pos_ls.append(r_start)
                        fixed_LTR_pos_ls.append(r_end)

                if len(fixed_LTR_bed) >= 2:
                    out_chrom = chrom
                    out_start = min(fixed_LTR_pos_ls)
                    out_end = max(fixed_LTR_pos_ls)
                
                elif len(fixed_LTR_bed) == 1:
                    out_chrom = chrom
                    for i in fixed_LTR_bed:
                        if (i + 1) in hap12_rmsk_LTR_dic[chrom] and hap12_rmsk_LTR_dic[chrom][(i + 1)][2] == fixed_LTR_bed[i][2] and hap12_rmsk_LTR_dic[chrom][(i + 1)][4] == fixed_LTR_bed[i][4]:
                            # print(hap12_rmsk_LTR_dic[chrom][i+1])
                            fixed_LTR_pos_ls.append(hap12_rmsk_LTR_dic[chrom][(i + 1)][0])
                            fixed_LTR_pos_ls.append(hap12_rmsk_LTR_dic[chrom][(i + 1)][1])
                            out_start = min(fixed_LTR_pos_ls)
                            out_end = max(fixed_LTR_pos_ls)
                        elif (i - 1) in hap12_rmsk_LTR_dic[chrom] and hap12_rmsk_LTR_dic[chrom][(i - 1)][2] == fixed_LTR_bed[i][2] and hap12_rmsk_LTR_dic[chrom][(i - 1)][4] == fixed_LTR_bed[i][4]:
                            fixed_LTR_pos_ls.append(hap12_rmsk_LTR_dic[chrom][(i - 1)][0])
                            fixed_LTR_pos_ls.append(hap12_rmsk_LTR_dic[chrom][(i - 1)][1])
                            out_start = min(fixed_LTR_pos_ls)
                            out_end = max(fixed_LTR_pos_ls)
                        else:
                            out_start = start
                            out_end = end
                else:
                    out_chrom = chrom
                    out_start = start
                    out_end = end

            else:
                out_chrom = chrom
                out_start = start
                out_end = end

            out_line = out_chrom + "\t" + str(out_start) + "\t" + str(out_end) + "\t" + ll[3] + "\n"
            outf.write(out_line)

if __name__ == "__main__":
    main()
