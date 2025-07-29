import sys
import pysam

def main():
    in_fixed_vcf = sys.argv[1] # AAC_AUC.ccs.deepconsensus.asm.split.fq.melTE.mm2sr.s.paf.combine.tab1.TEI.AAC.specific.vcf
    in_hap1_bam = sys.argv[2] # AAC.ccs.deepconsensus.fq.hap1.F4.s.bam
    in_hap2_bam = sys.argv[3] # AAC.ccs.deepconsensus.fq.hap2.F4.s.bam
    intag = sys.argv[4]

    out_bed = ".".join(in_fixed_vcf.split("/")[-1].split(".")[:-1]) + ".fixed.asm.bed"

    tei_rnames_dic = read_in_fixed_vcf_get_reads(in_fixed_vcf, intag)

    tei_pos_type_range_dic = {}
    chrom_count_dic = {}

    chrom_count_dic, tei_pos_type_range_dic = read_hap_bam(in_hap1_bam, tei_rnames_dic, chrom_count_dic, tei_pos_type_range_dic)
    chrom_count_dic, tei_pos_type_range_dic = read_hap_bam(in_hap2_bam, tei_rnames_dic, chrom_count_dic, tei_pos_type_range_dic)
    output_chrom_count_tei_pos_type_range(chrom_count_dic, tei_pos_type_range_dic, out_bed)


def read_in_fixed_vcf_get_reads(in_fixed_vcf, intag):
    tei_rnames_dic = {}
    with open(in_fixed_vcf, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")

            if ll[5] == "fixed":
                tei_pos_type = ll[0] + ":" + ll[1] + ":" + ll[2]

                if intag == "AAC" and ll[4] == "NA":
                    for r in ll[3].split(","):
                        tei_rnames_dic[r] = tei_pos_type

                if intag == "AUC" and ll[3] == "NA":
                    for r in ll[4].split(","):
                        tei_rnames_dic[r] = tei_pos_type

    return tei_rnames_dic

def read_hap_bam(in_bam, tei_rnames_dic, chrom_count_dic, tei_pos_type_range_dic):
    with pysam.AlignmentFile(in_bam, "rb") as samfile:
        for read in samfile:
            if read.query_name in tei_rnames_dic and not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                c = 0
                for cigar_type, cigar_length in read.cigartuples:
                    if (cigar_type == 4 or cigar_type == 5) and cigar_length >= 10 :  # Soft clip / Hard clip (clipped sequence present in SEQ)
                        c += 1
                        break

                if c == 0:
                    tei_pos_type = tei_rnames_dic[read.query_name]
                    chrom = read.reference_name
                    start = read.reference_start
                    end = read.reference_end

                    if tei_pos_type in chrom_count_dic and chrom in chrom_count_dic[tei_pos_type]:
                        chrom_count_dic.setdefault(tei_pos_type, {})[chrom] += 1
                    else:
                        chrom_count_dic.setdefault(tei_pos_type, {})[chrom] = 1

                    if tei_pos_type in tei_pos_type_range_dic and chrom in tei_pos_type_range_dic[tei_pos_type]:
                        tei_pos_type_range_dic[tei_pos_type][chrom] = (min(start, tei_pos_type_range_dic[tei_pos_type][chrom][0]), 
                                                                       max(end, tei_pos_type_range_dic[tei_pos_type][chrom][1]))
                    else:
                        tei_pos_type_range_dic.setdefault(tei_pos_type, {})[chrom] = (start, end)

    return chrom_count_dic, tei_pos_type_range_dic

def output_chrom_count_tei_pos_type_range(chrom_count_dic, tei_pos_type_range_dic, out_bed):
    with open(out_bed, "w") as outf:
        for tei_pos_type in chrom_count_dic:
            max_value = max(chrom_count_dic[tei_pos_type].values())
            for k, v in chrom_count_dic[tei_pos_type].items():
                if v == max_value:
                    most_depth_chrom = k
                    break
            
            if tei_pos_type in tei_pos_type_range_dic and most_depth_chrom in tei_pos_type_range_dic[tei_pos_type]:
                start, end = tei_pos_type_range_dic[tei_pos_type][most_depth_chrom]
                outline = most_depth_chrom + "\t" + str(start) + "\t" + str(end) + "\t" + tei_pos_type
                outf.write(outline + "\n")

if __name__ == "__main__":
    main()
