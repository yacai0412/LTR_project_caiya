import sys
import pysam
import gzip


def main():
    insnv = sys.argv[1] # AA.deepconsensus.cx3.fwd.retrogene.split.gene.fwd.AAC.mm2_splice.F904.shared.hap1.fwd_plus_minus.subpass_5.vaf_0.8.rm_25bp.snv
    dm6_bam = sys.argv[2] # AA.deepconsensus.cx3.fwd.retrogene.split.gene.dm6.mm2_splice.F904.s.bam
    fastq = sys.argv[3] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/TE_gene/retrogene/AA.deepconsensus.cx3.fwd.retrogene.split.gene.fastq

    out = insnv + ".consensuspos.out"

    snv_dic = read_bothstrand_fwd_snv(insnv)
    fastq_len_dic = read_fastq_get_rlen(fastq)
    read_bam_get_consnesus_pairwise_pos(dm6_bam, snv_dic, fastq_len_dic, out)


def read_bothstrand_fwd_snv(fwd_snv):
    snv_dic = {}
    with open(fwd_snv, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            support_read = ll[6]
            sep_name = support_read.split(":")[0]
            query_pos = int(support_read.split(":")[1])
            snv_dic.setdefault(sep_name, {})[query_pos] = line
    return snv_dic

def read_fastq_get_rlen(fastq):
    fastq_len_dic = {}
    c = 0
    
    # Check if the file is gzipped
    if fastq.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'  # Read as text
    else:
        open_func = open
        mode = 'r'   # Read as regular text file
    
    with open_func(fastq, mode) as inf:
        for line in inf:
            c += 1
            if c % 4 == 1:
                seqname = line.rstrip().lstrip("@")
            elif c % 4 == 2:
                rlen = len(line.rstrip())
                fastq_len_dic[seqname] = rlen
                
    return fastq_len_dic

def read_bam_get_consnesus_pairwise_pos(consensus_bam, snv_dic, fastq_len_dic, out):
    with pysam.AlignmentFile(consensus_bam, "rb") as samfile, open(out, "w") as final:
        for read in samfile:
            qname0 = read.query_name
            refname = read.reference_name
            if qname0 in snv_dic:
                fastq_length = fastq_len_dic[qname0]
                aligned_pairs_dic, strand = get_aligned_pairs_dic(read, fastq_length)
                for readpos in snv_dic[qname0]:
                    if readpos in aligned_pairs_dic:
                        consensus_pos = aligned_pairs_dic[readpos]
                    else:
                        closest_pos = min(aligned_pairs_dic.keys(), key=lambda k: abs(k - readpos))
                        consensus_pos = aligned_pairs_dic[closest_pos]

                    out_line = snv_dic[qname0][readpos] + "\t" + refname + "\t" + str(consensus_pos) + "\t" + strand
                    final.write(out_line + "\n")


def get_aligned_pairs_dic(read, fastq_length):
    aligned_pairs_tuple_ls = read.get_aligned_pairs(matches_only = True, with_seq = False)

    if read.cigartuples[0][0] == 5:
        left_hard_clip = read.cigartuples[0][1]
    else:
        left_hard_clip = 0

    if read.is_forward:
        aligned_pairs_dic = {key + left_hard_clip + 1: value + 1 for key, value in aligned_pairs_tuple_ls}
        strand = "+"
    elif read.is_reverse:
        aligned_pairs_dic = {fastq_length - (key + left_hard_clip + 1) + 1:  value + 1 for key, value in aligned_pairs_tuple_ls}
        strand = "-"
    else:
        print("mapping direction error in " + read.query_name)

    return aligned_pairs_dic, strand



if __name__ == "__main__":
    main()
