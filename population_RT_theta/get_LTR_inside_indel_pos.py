import sys
import pysam

def main():
    bam = sys.argv[1]

    out = ".".join(bam.split("/")[-1].split(".")[:-1]) + ".readpos.out"

    bam_ltr_inside_indel_pos_dic = read_bamfile(bam)
    write_inside_indel_pos(bam_ltr_inside_indel_pos_dic, out)


def read_bamfile(bam_file_path):
    bam_ltr_inside_indel_pos_dic = {}
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        for read in bam_file.fetch():
            rlen = read.query_length
            query_name = read.query_name
            target_name = read.reference_name            
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

                # Update read position based on CIGAR operation
                if op in [7, 8, 4]:  # Match Mismatch Softclip
                    read_pos += length
    
    return bam_ltr_inside_indel_pos_dic

def write_inside_indel_pos(bam_ltr_inside_indel_pos_dic, out):
    with open(out, "w") as final:
        for readname in bam_ltr_inside_indel_pos_dic:
            for indel_pos in bam_ltr_inside_indel_pos_dic[readname]:
                out_line = readname + "\t" + indel_pos
                final.write(out_line + "\n")

if __name__ == "__main__":
    main()

