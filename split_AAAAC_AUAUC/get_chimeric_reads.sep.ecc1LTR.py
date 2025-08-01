import sys
import gzip

def main():
    intab = sys.argv[1] # AA.deepconsensus.cx3.fwd.trans1_TE1.s.paf.te_gene.tab
    # no_split_readname = sys.argv[2] # /rd2/caiya/LTR/split_AAAAC_AUAUC/AA.deepconsensus.hap1hap2.nosplit.readnames
    infastqgz = sys.argv[2] # /rd/caiya/LTR/duplex/AA.deepconsensus.cx3.fwd.fastq.gz / /rd/caiya/LTR/duplex/AAC_ref/AA.deepconsensus.cx3.q30.fastq.gz

    fastq_out = ".".join(infastqgz.split("/")[-1].split(".")[:-2]) + ".ecc1LTR1.split.fastq"
    tab_out = ".".join(intab.split("/")[-1].split(".")[:-1]) + ".all.q30.tab"
    bed_out = ".".join(intab.split("/")[-1].split(".")[:-1]) + ".all.q30.bed"

    seq_dict = parse_fastq_gz(infastqgz)
    # no_split_readname_dic = get_nosplit_readnames_dic(no_split_readname)
    fwd_rev_paired_dic = read_intab_get_ecc_subfastq(intab)
    output_fwd_rev_paired_tab_fastq(fwd_rev_paired_dic, seq_dict, tab_out, fastq_out, bed_out)



def parse_fastq_gz(infastqgz):
    seq_dict = {}
    with gzip.open(infastqgz, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            f.readline()
            base_quality = f.readline().strip()
        
            seq_id = header.split("/")[0].lstrip("@") + "/" + header.split("/")[1] + "/deepconsensus/" + header.split("/")[3]            
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
    fwd_rev_paired_dic = {}
    with open(intab, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            rname = ll[0].split("/")[0] + "/" + ll[0].split("/")[1] + "/deepconsensus/" + ll[0].split("/")[3]
            zmw = ll[0].split("/")[1]
            query_strand = ll[0].split("/")[3]
            
            # if rname not in no_split_readname_dic:
            position_str = ll[2]
            cleaned_positions_ls, cleaned_positions_length_ls = remove_overlap(position_str)
            overlap_lengths_ls = calculate_overlap_lengths(position_str)
            if len(cleaned_positions_ls) > 1 and all(x >= 50 for x in cleaned_positions_length_ls):
                fwd_rev_paired_dic.setdefault(zmw, {})[query_strand] = rname, ll, cleaned_positions_ls, overlap_lengths_ls
    return fwd_rev_paired_dic

def output_fwd_rev_paired_tab_fastq(fwd_rev_paired_dic, seq_dict, tab_out, fastq_out, bed_out):
    with open(tab_out, "w") as tab_final, open(fastq_out, "w") as fastq_final, open(bed_out, "w") as bed_final:
        for zmw in fwd_rev_paired_dic:
            if len(fwd_rev_paired_dic[zmw]) == 2 and "fwd" in fwd_rev_paired_dic[zmw] and "rev" in fwd_rev_paired_dic[zmw]:
                for query_strand in fwd_rev_paired_dic[zmw]:
                    rname, ll, cleaned_positions_ls, overlap_lengths_ls = fwd_rev_paired_dic[zmw][query_strand]

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
                            fastq_final.write("@" + sub_name + "\n")
                            fastq_final.write(sub_seq + "\n")
                            fastq_final.write("+\n")
                            fastq_final.write(sub_baseq + "\n")

def remove_overlap(position_str):
    positions = []
    for pos in position_str.split(','):
        start, end = map(int, pos.split('-'))
        positions.append((start, end))

    cleaned_positions = []
    cleaned_positions_length = []
    prev_end = -1  # Initialize with a value before any possible start

    for start, end in positions:
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
