import sys
import pysam

def main():
    ref_fasta = sys.argv[1]
    inbam = sys.argv[2]
    query_fasta = sys.argv[3]

    snv_out = inbam.split("/")[-1] + ".snv"
    indel_out = inbam.split("/")[-1] + ".indel"

    ref_dic = read_fasta(ref_fasta)
    query_dic = read_fasta(query_fasta)
    read_bam(inbam, ref_dic, query_dic, snv_out, indel_out)


def read_fasta(ref_fasta):
    ref_dic = {}
    with open(ref_fasta, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line.lstrip(">")
                ref_dic[name] = ""
            else:
                ref_dic[name] = ref_dic[name] + line.upper()
    return ref_dic

def read_bam(inbam, ref_dic, query_dic, snv_out, indel_out):
    with pysam.AlignmentFile(inbam, "rb") as samfile, open(snv_out, "w") as final_snv, open(indel_out, "w") as final_indel:
        for read in samfile:
            # if not read.is_unmapped and not read.is_secondary:
            if not read.is_unmapped:
                mismatches, indels = get_mismatches_indels_from_cigar(read, query_dic, ref_dic)
                for alist in mismatches:
                    if alist[2] != "N":
                        final_snv.write("\t".join(alist) + "\n")
                
                for alist in indels:
                    if "N" not in alist[2]:
                        final_indel.write("\t".join(alist) + "\n")


def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    rev_complement = ''.join(complement[base] for base in reversed(dna_sequence))
    return rev_complement

# def get_mismatches_indels_from_cigar(read, fastq_len_dic, ref_dic):
def get_mismatches_indels_from_cigar(read, query_dic, ref_dic):
    mismatches = []
    indels = []
    query_pos = 0
    query_pos_total = 0
    ref_pos = read.reference_start
    qname0 = read.query_name

    if qname0 not in query_dic:
        print("error in fastq readlength " + qname0)
    else:
        rlen = len(query_dic[qname0])

        for cigar_type, cigar_length in read.cigartuples:
            if cigar_type == 7:  # Sequence match
                query_pos += cigar_length
                query_pos_total += cigar_length
                ref_pos += cigar_length
            elif cigar_type == 8:  # Sequence mismatch
                for i in range(cigar_length):
                    ref_base = ref_dic[read.reference_name][ref_pos]
                    # alt_base = read.query_sequence[query_pos]
                    alt_base = query_dic[qname0][query_pos_total]
                    query_pos += 1
                    query_pos_total += 1
                    ref_pos += 1

                    if read.is_forward:
                        query_pos1 = query_pos_total
                        rstrand = "+"
                    elif read.is_reverse:
                        query_pos1 = rlen - query_pos_total + 1
                        rstrand = "-"
                    qname = qname0 + ":" + str(query_pos1) + ":" + rstrand
                    mismatches.append((read.reference_name, str(ref_pos), ref_base, alt_base, qname))
                    
            elif cigar_type == 1:  # Insertion to the reference (query pos == including the insertion base previous base, pos same as the alt base)
                insertion_start = query_pos_total + 1
                ins_end = query_pos_total + cigar_length 
                ref_base = ref_dic[read.reference_name][ref_pos - 1]
                # ref_base = "-"
                # alt_base = read.query_sequence[query_pos - 1: query_pos + cigar_length]
                alt_base = query_dic[qname0][query_pos_total - 1: query_pos_total + cigar_length]

                if read.is_forward:
                    # query_pos1 = insertion_start
                    query_pos1 = query_pos_total
                    rstrand = "+"
                elif read.is_reverse:
                    # query_pos1 = rlen - ins_end + 1
                    query_pos1 = rlen - ins_end
                    rstrand = "-"
                qname = qname0 + ":" + str(query_pos1) + ":" + rstrand
                indels.append((read.reference_name, str(ref_pos), ref_base, alt_base, qname, str(cigar_length), "INS"))

                query_pos += cigar_length
                query_pos_total += cigar_length

            elif cigar_type == 2:  # Deletion from the reference (query pos == including the boundary 1bp query base)
                # h1tg000001l\t879952\tCT\tC\tm64054_211209_114633/19007732/deepconsensus/fwd_3:127:+\t1\tDEL, read 127 base is C, same as previous alt base
                ref_base = ref_dic[read.reference_name][ref_pos - 1: ref_pos + cigar_length]
                # alt_base = "-"
                # alt_base = read.query_sequence[query_pos - 1]
                alt_base = query_dic[qname0][query_pos_total - 1]

                if read.is_forward:
                    query_pos1 = query_pos_total
                    rstrand = "+"
                elif read.is_reverse:
                    query_pos1 = rlen - query_pos_total
                    rstrand = "-"

                qname = qname0 + ":" + str(query_pos1) + ":" + rstrand
                indels.append((read.reference_name, str(ref_pos), ref_base, alt_base, qname, str(cigar_length), "DEL"))
                ref_pos += cigar_length

            elif cigar_type == 3:  # BAM CREF SKIP (N) in minimap2/winnowmap2 splice mode
                ref_pos += cigar_length
            elif cigar_type == 4 :  # Soft clip (clipped sequence present in SEQ)
                query_pos += cigar_length
                query_pos_total += cigar_length
            elif cigar_type == 5:  # Hard clip (clipped sequence absent from SEQ)
                query_pos_total += cigar_length

    return mismatches, indels

# def get_mismatches_from_cigar(read, ref_dic, query_dic):
#     mismatches = []
#     query_pos = 0
#     ref_pos = read.reference_start
#     rlen = len(query_dic[read.query_name])
#     rq = read.get_tag("rq")

#     for cigar_type, cigar_length in read.cigartuples:
#         if cigar_type == 7:  # Sequence match
#             query_pos += cigar_length
#             ref_pos += cigar_length
#         elif cigar_type == 8:  # Sequence mismatch
#             for i in range(cigar_length):
#                 ref_base = ref_dic[read.reference_name][ref_pos]
#                 alt_base = query_dic[read.query_name][query_pos]
#                 query_pos += 1
#                 ref_pos += 1

#                 if read.is_forward:
#                     query_pos1 = query_pos
#                     rstrand = "+"
#                     alt_base = query_dic[read.query_name][query_pos1 - 1]
#                 elif read.is_reverse:
#                     query_pos1 = rlen - query_pos + 1
#                     rstrand = "-"
#                     alt_base = reverse_complement(query_dic[read.query_name][query_pos1 - 1])


#                 qname = read.query_name + ":" + str(query_pos1) + ":" + str(rq) + ":" + rstrand
#                 qname = read.query_name + ":" + str(query_pos1) + ":" + rstrand
#                 mismatches.append((read.reference_name, str(ref_pos), ref_base, alt_base, qname))
                
#         elif cigar_type == 1:  # Insertion to the reference
#             query_pos += cigar_length
#         elif cigar_type == 2:  # Deletion from the reference
#             ref_pos += cigar_length
#         elif cigar_type == 4:  # Soft clip (clipped sequence present in SEQ)
#             query_pos += cigar_length
#         elif cigar_type == 5:  # Hard clip (clipped sequence absent from SEQ)
#             pass

#     return mismatches



if __name__ == "__main__":
    main()


