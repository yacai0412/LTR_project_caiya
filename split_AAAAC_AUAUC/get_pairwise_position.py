import sys
import gzip
import pysam
import re

def main():
    fwd_snv = sys.argv[1]
    consensus_bam = sys.argv[2]
    fastq = sys.argv[3]
    ref = sys.argv[4]

    depth = consensus_bam + ".depth"
    out = fwd_snv + ".consensuspos.out"
    out1 = fwd_snv + ".consensuspos.out1"

    snv_dic = read_bothstrand_fwd_snv(fwd_snv)
    fastq_len_dic = read_fastq_get_rlen(fastq)
    ref_pos_dic = create_ref_length_pos_dic(ref)
    read_bam_get_consnesus_pairwise_pos(consensus_bam, snv_dic, fastq_len_dic, ref_pos_dic, out)

    depth_dic = read_depth(depth)
    export_ref_pos_count(ref_pos_dic, depth_dic, out1)

def create_ref_length_pos_dic(ref):
    ref_len_dic = {}
    with open(ref, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            if line.startswith(">"):
                rname = line.split("\t")[0].lstrip(">")
                ref_len_dic[rname] = 0
            else:
                ref_len_dic[rname] += len(line)

    ref_pos_dic = {}
    for rname in ref_len_dic:
        for i in range(ref_len_dic[rname]):
            ref_pos_dic.setdefault(rname, {})[i + 1] = 0
    return ref_pos_dic


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

def read_depth(depth):
    depth_dic = {}
    with open(depth, "r") as inf:
        for line in inf:
            ll = line.rstrip("\n").split("\t")
            depth_dic.setdefault(ll[0], {})[int(ll[1])] = int(ll[2])
    return depth_dic

def get_aligned_pairs_dic(read, fastq_length):
    aligned_pairs_tuple_ls = read.get_aligned_pairs(matches_only = True, with_seq = False)

    if read.cigartuples[0][0] == 5:
        left_hard_clip = read.cigartuples[0][1]
    else:
        left_hard_clip = 0

    if read.is_forward:
        aligned_pairs_dic = {key + left_hard_clip + 1: value + 1 for key, value in aligned_pairs_tuple_ls}
    elif read.is_reverse:
        aligned_pairs_dic = {fastq_length - (key + left_hard_clip + 1) + 1:  value + 1 for key, value in aligned_pairs_tuple_ls}
    else:
        print("mapping direction error in " + read.query_name)

    return aligned_pairs_dic


def read_bam_get_consnesus_pairwise_pos(consensus_bam, snv_dic, fastq_len_dic, ref_pos_dic, out):
    with pysam.AlignmentFile(consensus_bam, "rb") as samfile, open(out, "w") as final:
        for read in samfile:
            qname0 = read.query_name
            refname = read.reference_name
            if qname0 in snv_dic:
                fastq_length = fastq_len_dic[qname0]
                aligned_pairs_dic = get_aligned_pairs_dic(read, fastq_length)
                for readpos in snv_dic[qname0]:
                    if readpos in aligned_pairs_dic:
                        consensus_pos = aligned_pairs_dic[readpos]
                    else:
                        closest_pos = min(aligned_pairs_dic.keys(), key=lambda k: abs(k - readpos))
                        consensus_pos = aligned_pairs_dic[closest_pos]

                    ref_pos_dic[refname][consensus_pos] += 1
                    out_line = snv_dic[qname0][readpos] + "\t" + refname + "\t" + str(consensus_pos)
                    final.write(out_line + "\n")
    return ref_pos_dic


def export_ref_pos_count(ref_pos_dic, depth_dic, out1):
    with open(out1, "w") as final1:
        for ltr in ref_pos_dic:
            for ref_pos in ref_pos_dic[ltr]:
                count = ref_pos_dic[ltr][ref_pos]
                if ltr in depth_dic and ref_pos in depth_dic[ltr]:
                    depth = depth_dic[ltr][ref_pos]
                else:
                    depth = 0

                if depth != 0:
                    mutation_rate = count / depth
                else:
                    mutation_rate = 0

                if ltr == "DMLTR5":
                    ref_ltr_name = "HMSBEAGLE"
                    ref_ltr_type = "LTR"
                else:        
                    ref_ltr_name = "_".join(re.split(r'[-_]', ltr)[:-1])
                    ref_ltr_type = re.split(r'[-_]', ltr)[-1]

                out_line = ltr + "\t" + ref_ltr_name + "\t" + ref_ltr_type + "\t" + str(ref_pos) + "\t" + str(count) + "\t" + str(depth) + "\t" + str(mutation_rate)
                final1.write(out_line + "\n")



if __name__ == "__main__":
    main()
