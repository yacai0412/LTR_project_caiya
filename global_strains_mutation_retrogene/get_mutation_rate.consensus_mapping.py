import sys

def main():
    in_consensus_snv = sys.argv[1] # G0.G0-F100.G73.TEIseq.ont.wtdbg2.ctg.LTR_I_LTR_consensus.mm2ont.F904.s.bam.indel.G73.consensus_pos.count
    in_depth = sys.argv[2] # G0.G0-F100.G73.TEIseq.ont.wtdbg2.ctg.LTR_I_LTR_consensus.mm2ont.F904.s.G73.bam.depth

    out = in_consensus_snv.replace(".count", ".mutation.rate")

    depth_dic = read_in_depth(in_depth)
    calculate_mutation_rate(in_consensus_snv, depth_dic, out)

def read_in_depth(in_depth):
    depth_dic = {}
    with open(in_depth, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")

            chrom = ll[0]
            pos = int(ll[1])
            depth = int(ll[2])

            depth_dic.setdefault(chrom, {})[pos] = depth
    return depth_dic

def calculate_mutation_rate(in_consensus_snv, depth_dic, out):
    with open(in_consensus_snv, "r") as inf, open(out, "w") as outf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")

            chrom = ll[0]
            pos = int(ll[1])
            mutation_count = int(ll[2])

            # Get the depth for this position
            depth = depth_dic.get(chrom, {}).get(pos, 0)

            if depth > 0:
                mutation_rate = mutation_count / depth
            else:
                mutation_rate = 0.0

            outf.write(f"{chrom}\t{pos}\t{mutation_count}\t{depth}\t{mutation_rate}\n")


if __name__ == "__main__":
    main()
