import sys
import pysam


def main():
    hap1_bam = sys.argv[1]
    hap2_bam = sys.argv[2]

    out = sys.argv[3]

    hap1_no_split_dic = read_bam_get_no_split_reads(hap1_bam)
    hap2_no_split_dic = read_bam_get_no_split_reads(hap2_bam)
    out_put_no_split(hap1_no_split_dic, hap2_no_split_dic, out)


def read_bam_get_no_split_reads(inbam):
    no_split_dic = {}
    with pysam.AlignmentFile(inbam, "rb") as samfile:
        for read in samfile:
            if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                rname = read.query_name
                cigar_tuples = read.cigartuples
                ciagr_string = read.cigarstring
                if (cigar_tuples[0][0] != 5 and cigar_tuples[0][0] != 4) and (cigar_tuples[-1][0] != 5 and cigar_tuples[-1][0] != 4):
                    no_split_dic[rname] = ciagr_string
    return no_split_dic

def out_put_no_split(hap1_no_split_dic, hap2_no_split_dic, out):
    total_no_split_dic = {}
    for rname1 in hap1_no_split_dic:
        if rname1 in hap2_no_split_dic:
            total_no_split_dic[rname1] = "hap1_hap2\t" + hap1_no_split_dic[rname1] + " " + hap2_no_split_dic[rname1]
        else:
            total_no_split_dic[rname1] = "hap1\t" + hap1_no_split_dic[rname1]
    
    for rname2 in hap2_no_split_dic:
        if rname2 not in total_no_split_dic:
            total_no_split_dic[rname2] = "hap2\t" + hap2_no_split_dic[rname2]
    
    with open(out, "w") as final:
        for rname0 in total_no_split_dic:
            outline = rname0 + "\t" + total_no_split_dic[rname0]
            final.write(outline + "\n")


if __name__ == "__main__":
    main()
