import sys
import pysam
import gzip

def main():
    ref_fasta = sys.argv[1] # /rd/caiya/dm6.major.chr.fasta
    inbam = sys.argv[2]
    fastq = sys.argv[3]
    out = inbam.split("/")[-1] + ".snv"
    final = open(out, "w")

    ref_dic = read_fasta(ref_fasta)
    fastq_len_dic = read_fastq_get_rlen(fastq)
    read_bam(inbam, fastq_len_dic, ref_dic, final)
    final.close()


def read_fasta(ref_fasta):
    ref_dic = {}
    with open(ref_fasta, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line.lstrip(">")
                ref_dic[name] = ""
            else:
                ref_dic[name] = ref_dic[name] + line    
    return ref_dic

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

def read_bam(inbam, fastq_len_dic, ref_dic, final):
    with pysam.AlignmentFile(inbam, "rb") as samfile:
        for read in samfile:
            if not read.is_unmapped and not read.is_secondary:
                mismatches = get_mismatches_from_cigar(read, fastq_len_dic, ref_dic)
                for alist in mismatches:
                    if alist[2] != "N":
                        final.write("\t".join(alist) + "\n")

def get_mismatches_from_cigar(read, fastq_len_dic, ref_dic):
    mismatches = []
    query_pos = 0
    query_pos_total = 0
    ref_pos = read.reference_start
    qname0 = read.query_name

    if qname0 not in fastq_len_dic:
        print("error in fastq readlength " + qname0)
    else:
        rlen = fastq_len_dic[qname0]

        for cigar_type, cigar_length in read.cigartuples:
            if cigar_type == 7:  # Sequence match
                query_pos += cigar_length
                query_pos_total += cigar_length
                ref_pos += cigar_length
            elif cigar_type == 8:  # Sequence mismatch
                for i in range(cigar_length):
                    ref_base = ref_dic[read.reference_name][ref_pos]
                    alt_base = read.query_sequence[query_pos]
                    query_pos += 1
                    query_pos_total += 1
                    ref_pos += 1

                    if read.is_forward:
                        query_pos1 = query_pos_total
                        rstrand = "+"
                    elif read.is_reverse:
                        query_pos1 = rlen - query_pos_total + 1
                        rstrand = "-"
                    # qname = read.query_name + ":" + str(query_pos1) + ":" + str(rq) + ":" + rstrand
                    qname = read.query_name + ":" + str(query_pos1) + ":" + rstrand
                    mismatches.append((read.reference_name, str(ref_pos), ref_base, alt_base, qname))
                    
            elif cigar_type == 1:  # Insertion to the reference
                query_pos += cigar_length
                query_pos_total += cigar_length
            elif cigar_type == 2 or cigar_type == 3:  # Deletion from the reference / BAM CREF SKIP (N) in minimap2/winnowmap2 splice mode
                ref_pos += cigar_length
            elif cigar_type == 4 :  # Soft clip (clipped sequence present in SEQ)
                query_pos += cigar_length
                query_pos_total += cigar_length
            elif cigar_type == 5:  # Hard clip (clipped sequence absent from SEQ)
                query_pos_total += cigar_length

    return mismatches

if __name__ == "__main__":
    main()


