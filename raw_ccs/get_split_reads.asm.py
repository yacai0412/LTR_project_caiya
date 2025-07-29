import sys
import gzip

def main():
    intab = sys.argv[1] # AA.ori.ccs.fq.AAC.hap1.win.s.uniq.paf.asm.combine.shared.chimeric.tab
    infastqgz = sys.argv[2] # AA.ori.ccs.fastq.gz
    out_tag = sys.argv[3] # asm

    if infastqgz.endswith(".fastq.gz"):
        fastq_out = ".".join(infastqgz.split("/")[-1].split(".")[:-2]) + "." + out_tag + ".split.fastq"
    elif infastqgz.endswith(".fastq"):
        fastq_out = ".".join(infastqgz.split("/")[-1].split(".")[:-1]) + "." + out_tag + ".split.fastq"

    tab_out = ".".join(intab.split("/")[-1].split(".")[:-1]) + "." + out_tag + ".split.tab"
    bed_out = ".".join(intab.split("/")[-1].split(".")[:-1]) + "." + out_tag + ".split.bed"

    seq_dict = parse_fastq_gz(infastqgz)
    oriccs_dic = read_intab_get_ecc_subfastq(intab)
    output_oriccs_tab_fastq(oriccs_dic, seq_dict, tab_out, fastq_out, bed_out)


def parse_fastq_gz(infastqgz):
    seq_dict = {}

    if infastqgz.endswith(".fastq"):
        open_func = open
        mode = "r"
    elif infastqgz.endswith(".fastq.gz"):
        open_func = gzip.open
        mode = "rt"

    with open_func(infastqgz, mode) as f:
        while True:
            header = f.readline().strip().split(" ")[0]
            if not header:
                break
            sequence = f.readline().strip()
            f.readline()
            base_quality = f.readline().strip()
        
            seq_id = header.lstrip("@")            
            seq_dict[seq_id] = (sequence, base_quality)
    return seq_dict

def get_nosplit_readnames_dic(no_split_readname):
    no_split_readname_dic = {}
    with open(no_split_readname, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            rname = ll[0]
            no_split_readname_dic[rname] = 1
    return no_split_readname_dic

def read_intab_get_ecc_subfastq(intab):
    oriccs_dic = {}
    with open(intab, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            rname = ll[0]
            
            # if rname not in no_split_readname_dic:
            position_str = ll[2]
            cleaned_positions_ls, cleaned_positions_length_ls = remove_overlap(position_str)
            overlap_lengths_ls = calculate_overlap_lengths(position_str)
            if len(cleaned_positions_ls) > 1 and all(x >= 50 for x in cleaned_positions_length_ls):
                oriccs_dic[rname] = ll, cleaned_positions_ls, overlap_lengths_ls
    return oriccs_dic

def output_oriccs_tab_fastq(oriccs_dic, seq_dict, tab_out, fastq_out, bed_out):
    with open(tab_out, "w") as tab_final, open(fastq_out, "w") as fastq_final, open(bed_out, "w") as bed_final:
        for rname in oriccs_dic:
            ll, cleaned_positions_ls, overlap_lengths_ls = oriccs_dic[rname]

            if rname in seq_dict:
                tab_final.write(rname + "\t" + "\t".join(ll[1:]) + "\t" + ",".join(overlap_lengths_ls) + "\n")
                c = 0
                for pos_tuples in cleaned_positions_ls:
                    c += 1
                    start = int(pos_tuples[0])
                    end = int(pos_tuples[1])

                    sub_name = rname + "_" + str(c)
                    bed_out_line = rname + "\t" + str(start) + "\t" + str(end) + "\t" + sub_name
                    bed_final.write(bed_out_line + "\n")

                    sub_seq = seq_dict[rname][0][start: end]
                    sub_baseq = seq_dict[rname][1][start: end]
                    
                    out_sub_fastq = "@" + sub_name + "\n" + sub_seq + "\n+\n" + sub_baseq + "\n"
                    fastq_final.write(out_sub_fastq)

# def remove_full_overlap(position_str):
#     positions = []
#     for pos in position_str.split(','):
#         start, end = map(int, pos.split('-'))
#         positions.append((start, end))
#     positions.sort(key=lambda x: (x[0], -x[1]))
#     filtered_positions = []
#     max_end = -1
#     for start, end in positions:
#         if end > max_end:
#             max_end = end
#             filtered_positions.append((start, end))
#         else:
#             continue
#     return filtered_positions

def remove_overlap(position_str):
    # filtered_positions = remove_full_overlap(position_str)
    position_ls = position_str.split(',')
    
    cleaned_positions = []
    cleaned_positions_length = []
    prev_end = -1  # Initialize with a value before any possible start

    for pos0 in position_ls:
        start = int(pos0.split("-")[0])
        end = int(pos0.split("-")[1])
        if start < prev_end:
            start = prev_end

        cleaned_positions.append((start, end))
        cleaned_positions_length.append(end - start)
        prev_end = end

    return cleaned_positions, cleaned_positions_length

def calculate_overlap_lengths(position_str):
    positions = []
    for pos in position_str.split(','):
        start, end = map(int, pos.split('-'))
        positions.append((start, end))

    overlap_lengths = []

    for i in range(1, len(positions)):
        prev_start, prev_end = positions[i - 1]
        curr_start, curr_end = positions[i]

        if curr_start < prev_end:
            overlap_length = min(prev_end, curr_end) - curr_start
            overlap_lengths.append(str(overlap_length))
        else:
            overlap_lengths.append(str(0))  # No overlap

    return overlap_lengths

if __name__ == "__main__":
    main()
