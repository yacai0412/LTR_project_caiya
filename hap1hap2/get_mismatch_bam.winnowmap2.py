import sys
import pysam

def main():
    ref_fasta = sys.argv[1] # /rd/caiya/dm6.major.chr.fasta
    inbam = sys.argv[2]
    out = inbam.split("/")[-1] + ".snv"
    final = open(out, "w")

    ref_dic = read_fasta(ref_fasta)
    read_bam(inbam, ref_dic, final)
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

def read_bam(inbam, ref_dic, final):
    with pysam.AlignmentFile(inbam, "rb") as samfile:
        for read in samfile:
            if not read.is_unmapped and not read.is_secondary:
                mismatches = get_mismatches_from_cigar(read, ref_dic)
                for alist in mismatches:
                    if alist[2] != "N":
                        final.write("\t".join(alist) + "\n")

def get_mismatches_from_cigar(read, ref_dic):
    mismatches = []
    query_pos = 0
    ref_pos = read.reference_start
    rlen = read.query_length
    # rq = read.get_tag("rq")

    for cigar_type, cigar_length in read.cigartuples:
        if cigar_type == 7:  # Sequence match
            query_pos += cigar_length
            ref_pos += cigar_length
        elif cigar_type == 8:  # Sequence mismatch
            for i in range(cigar_length):
                ref_base = ref_dic[read.reference_name][ref_pos]
                alt_base = read.query_sequence[query_pos]
                query_pos += 1
                ref_pos += 1

                if read.is_forward:
                    query_pos1 = query_pos
                    rstrand = "+"
                elif read.is_reverse:
                    query_pos1 = rlen - query_pos + 1
                    rstrand = "-"
                # qname = read.query_name + ":" + str(query_pos1) + ":" + str(rq) + ":" + rstrand
                qname = read.query_name + ":" + str(query_pos1) + ":" + rstrand
                mismatches.append((read.reference_name, str(ref_pos), ref_base, alt_base, qname))
                
        elif cigar_type == 1:  # Insertion to the reference
            query_pos += cigar_length
        elif cigar_type == 2:  # Deletion from the reference
            ref_pos += cigar_length
        elif cigar_type == 4:  # Soft clip (clipped sequence present in SEQ)
            query_pos += cigar_length
        elif cigar_type == 5:  # Hard clip (clipped sequence absent from SEQ)
            pass

    return mismatches

if __name__ == "__main__":
    main()


