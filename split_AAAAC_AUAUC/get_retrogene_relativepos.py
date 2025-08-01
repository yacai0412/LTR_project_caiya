import sys

def main():
    dm6_largest_transcript = "/rd/caiya/dm6.refGene.largest_transcript.gtf"

    dm6_consensus_position = sys.argv[1]
    out = ".".join(dm6_consensus_position.split("/")[-1].split(".")[:-1]) + ".gene_relative_pos.out"
    out1 = ".".join(dm6_consensus_position.split("/")[-1].split(".")[:-1]) + ".gene_relative_pos.out1"


    dm6_transcrpt_pos_dic = read_dm6_largest_transcript(dm6_largest_transcript)
    relative_pos1_ls = read_dm6_consensus_position_get_gene_pos(dm6_consensus_position, dm6_transcrpt_pos_dic, out)
    get_realtive_pos1_ref(relative_pos1_ls, out1)


def read_dm6_largest_transcript(dm6_largest_transcript):
    dm6_transcrpt_pos_dic = {}
    with open(dm6_largest_transcript, "r") as inf:
        for line in inf:
            ll = line.rstrip("\n").split("\t")
            chrom = ll[0]
            start_end = ll[3] + "\t" + ll[4]
            strand = ll[6]
            gene_id = ll[8].split(";")[0].split(" ")[1].strip('"')

            dm6_transcrpt_pos_dic.setdefault(chrom, {})[start_end] = strand, gene_id
    return dm6_transcrpt_pos_dic

def read_dm6_consensus_position_get_gene_pos(dm6_consensus_position, dm6_transcrpt_pos_dic, out):
    relative_pos1_ls = []
    with open(dm6_consensus_position, "r") as inf, open(out, "w") as final:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")

            dm6_chrom = ll[-3]
            dm6_pos = int(ll[-2])

            if dm6_chrom in dm6_transcrpt_pos_dic:
                relative_gene_id = "NA"
                relative_gene_pos = "NA"
                relative_gene_pos1 = "NA"
                for start_end in dm6_transcrpt_pos_dic[dm6_chrom]:
                    start = int(start_end.split("\t")[0])
                    end = int(start_end.split("\t")[1])

                    if start <= dm6_pos <= end:
                        gene_strand, gene_id = dm6_transcrpt_pos_dic[dm6_chrom][start_end]
                        gene_length = end - start + 1

                        relative_gene_id = gene_id
                        if gene_strand == "+":
                            relative_gene_pos = dm6_pos - start + 1
                        elif gene_strand == "-":
                            relative_gene_pos = end - dm6_pos + 1
                        else:
                            print("gene strand error in " + gene_strand + "\t" + gene_id)

                        if relative_gene_pos <= gene_length / 2:
                            relative_gene_pos1 = relative_gene_pos
                        else:
                            relative_gene_pos1 = -1 * (gene_length - relative_gene_pos + 1)

                        break
                
                relative_pos1_ls.append(relative_gene_pos1)
                out_line = line + "\t" + relative_gene_id + "\t" + str(relative_gene_pos) + "\t" + str(relative_gene_pos1)
                final.write(out_line + "\n")

    return relative_pos1_ls

def get_realtive_pos1_ref(relative_pos1_ls, out1):
    relative_pos1_ref_dic = {}
    max_relativepos1 = max(relative_pos1_ls)
    min_relativepos1 = min(relative_pos1_ls)

    max_abs_relativepos1 = max(abs(max_relativepos1), abs(min_relativepos1))
    relative_pos1_ref_dic[0] = 0
    for i in range(max_abs_relativepos1):
        relative_pos1_ref_dic[i + 1] = 0
        relative_pos1_ref_dic[-1 * (i + 1)] = 0

    for pos1 in relative_pos1_ls:
        relative_pos1_ref_dic[pos1] += 1
    
    with open(out1, "w") as final1:
        for pos11 in relative_pos1_ref_dic:
            out_line = str(pos11) + "\t" + str(relative_pos1_ref_dic[pos11])
            final1.write(out_line + "\n")



if __name__ == "__main__":
    main()
