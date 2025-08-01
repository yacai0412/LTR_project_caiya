import sys
import gzip

def main():
    insnv = sys.argv[1] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/both_strand/no_filter_softclip/denovo_mutation/AA.indel.consensuspos.out
    infastq = sys.argv[2] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/ecc/AA.deepconsensus.cx3.q30.ecc1LTR1.split.fastq
    small_reads_boundary = int(sys.argv[3]) # 25

    out_fastq = ".".join(insnv.split("/")[-1].split(".")[:-1]) + "." + str(small_reads_boundary) + "bp.fastq"

    snv_query_pos_dic = read_snv(insnv, small_reads_boundary)
    fastq_seq_dic = read_fastq_get_seq(infastq)

    get_small_reads(snv_query_pos_dic, fastq_seq_dic, out_fastq)

def read_snv(insnv, small_reads_boundary):
    snv_query_pos_dic = {}
    with open(insnv, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            snv_id = ll[0]
            ref_base = ll[3]
            alt_base = ll[4]

            query_name = ll[5].split(":")[0]
            query_pos = int(ll[5].split(":")[1])
            query_start = query_pos - small_reads_boundary
            query_end = query_pos + max(0, (len(alt_base) - len(ref_base))) + small_reads_boundary
            snv_query_pos_dic[snv_id] = query_name, query_start, query_end
    return snv_query_pos_dic


def read_fastq_get_seq(fastq):
    fastq_seq_dic = {}
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
            line = line.rstrip("\n")
            c += 1
            if c % 4 == 1:
                seqname = line.lstrip("@")
            elif c % 4 == 2:
                fastq_seq_dic[seqname] = line
            elif c % 4 == 0:
                fastq_seq_dic[seqname] = fastq_seq_dic[seqname] + "\t" + line

    return fastq_seq_dic


def get_small_reads(snv_query_pos_dic, fastq_seq_dic, out_fastq):
    with open(out_fastq, "w") as final_fastq:
        for snv_id in snv_query_pos_dic:
            query_name, query_start, query_end = snv_query_pos_dic[snv_id]
            fastq_seq_quality = fastq_seq_dic[query_name]
            fastq_seq = fastq_seq_quality.split("\t")[0]
            fastq_quality = fastq_seq_quality.split("\t")[1]
            len_fastq_seq = len(fastq_seq)

            subseq = fastq_seq[max(query_start - 1, 0) : min(query_end, len_fastq_seq)]
            subquality = fastq_quality[max(query_start - 1, 0) : min(query_end, len_fastq_seq)]

            out_line = "@" + snv_id + "\n" + subseq + "\n+\n" + subquality + "\n"
            final_fastq.write(out_line)


if __name__ == "__main__":
    main()
