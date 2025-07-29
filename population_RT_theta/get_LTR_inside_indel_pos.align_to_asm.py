import sys
import pysam

def main():
    bam = sys.argv[1]
    infasta = sys.argv[2] # 

    out = ".".join(bam.split("/")[-1].split(".")[:-1]) + ".readpos.out"
    out1 = ".".join(bam.split("/")[-1].split(".")[:-1]) + ".readpos.specific.out"
    out2 = ".".join(bam.split("/")[-1].split(".")[:-1]) + ".readpos.specific.largest.out"

    fasta_dic = read_fasta(infasta)
    bam_ltr_inside_indel_pos_dic, bam_ltr_inside_indel_pos_count_dic = read_bamfile(bam, fasta_dic)
    write_inside_indel_pos(bam_ltr_inside_indel_pos_dic, bam_ltr_inside_indel_pos_count_dic, out, out1, out2)


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

def read_bamfile(bam_file_path, fasta_dic):
    bam_ltr_inside_indel_pos_dic = {}
    bam_ltr_inside_indel_pos_count_dic = {}
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        for read in bam_file.fetch():
            # rlen = read.query_length
            query_name = read.query_name
            target_name = read.reference_name            
            rlen = len(fasta_dic[query_name])
            query_target_name = query_name + "\t" + target_name
            read_pos = 0

            for op, length in read.cigartuples:
                if op == 1: # insertion
                    insertion_start = read_pos + 1
                    ins_end = read_pos + length 

                    if read.is_forward:
                        query_pos1 = insertion_start
                        rstrand = "+"
                    elif read.is_reverse:
                        query_pos1 = rlen - ins_end + 1
                        rstrand = "-"

                    indel_pos = str(query_pos1) + ":" + rstrand + ":" + str(length) + "bp_insertion"
                    bam_ltr_inside_indel_pos_dic.setdefault(query_target_name, []).append(indel_pos)
                    if query_name in bam_ltr_inside_indel_pos_count_dic and indel_pos in bam_ltr_inside_indel_pos_count_dic[query_name]:
                        bam_ltr_inside_indel_pos_count_dic[query_name][indel_pos] += 1
                    else:
                        bam_ltr_inside_indel_pos_count_dic.setdefault(query_name, {})[indel_pos] = 1
                    read_pos = ins_end

                elif op == 2: # deletion
                    if read.is_forward:
                        query_pos1 = read_pos
                        rstrand = "+"
                    elif read.is_reverse:
                        query_pos1 = rlen - read_pos
                        rstrand = "-"

                    indel_pos = str(query_pos1) + ":" + rstrand + ":" + str(length) + "bp_deletion"
                    bam_ltr_inside_indel_pos_dic.setdefault(query_target_name, []).append(indel_pos)
                    if query_name in bam_ltr_inside_indel_pos_count_dic and indel_pos in bam_ltr_inside_indel_pos_count_dic[query_name]:
                        bam_ltr_inside_indel_pos_count_dic[query_name][indel_pos] += 1
                    else:
                        bam_ltr_inside_indel_pos_count_dic.setdefault(query_name, {})[indel_pos] = 1

                # Update read position based on CIGAR operation
                if op in [7, 8, 4]:  # Match Mismatch Softclip
                    read_pos += length
    
    return bam_ltr_inside_indel_pos_dic, bam_ltr_inside_indel_pos_count_dic

def write_inside_indel_pos(bam_ltr_inside_indel_pos_dic, bam_ltr_inside_indel_pos_count_dic, out, out1, out2):
    with open(out, "w") as final, open(out1, "w") as final1, open(out2, "w") as final2:
        for readname in bam_ltr_inside_indel_pos_dic:
            for indel_pos in bam_ltr_inside_indel_pos_dic[readname]:
                out_line = readname + "\t" + indel_pos
                final.write(out_line + "\n")

        for query_name in bam_ltr_inside_indel_pos_count_dic:
            out_query_name = query_name.split(":")[0] + "\t" + query_name.split(":")[1] + ":" + query_name.split(":")[2]
            max_length = 0
            max_length_indel_pos = ""
            for indel_pos in bam_ltr_inside_indel_pos_count_dic[query_name]:
                if bam_ltr_inside_indel_pos_count_dic[query_name][indel_pos] == 1:
                    indel_length = int(indel_pos.split(":")[2].split("bp_")[0])
                    if indel_length >= max_length:
                        max_length = indel_length
                        max_length_indel_pos = indel_pos
                    
                    out_line1 = out_query_name + "\t" + indel_pos
                    final1.write(out_line1 + "\n")

            if max_length_indel_pos != "":
                out_line2 = out_query_name + "\t" + max_length_indel_pos
                final2.write(out_line2 + "\n")



if __name__ == "__main__":
    main()

