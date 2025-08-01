import sys
import re

def main():
    consensuspos_snv = sys.argv[1] # AA.deepconsensus.cx3.q30.ecc1LTR1.no_filter_softclip.split.fwd.AAC.F904.shared.hap1.fwd_plus_minus.subpass_5.vaf_0.8.rm_25bp.consensuspos.out / AA.deepconsensus.cx3.q30.ecc1LTR1.no_filter_softclip.split.fwd.AAC.F904.shared.hap1.indel.fwd_plus_minus.subpass_8.vaf_0.8.rm_25bp.1bp_DEL.consensuspos.out
    depth = sys.argv[2] # mapping to consensus
    ref = sys.argv[3]

    out1 = ".".join(consensuspos_snv.split(".")[:-1]) + ".mutation_type.out1"

    snv_count_dic, ref_consensus_pos_dic = read_consensuspos_snv(consensuspos_snv)
    consensus_pos_count_dic = create_ref_length_pos_dic(ref)
    depth_dic = read_depth(depth)
    export_ref_pos_count_mutation_type(snv_count_dic, ref_consensus_pos_dic, consensus_pos_count_dic, depth_dic, out1)


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

    consensus_pos_count_dic = {}
    for rname in ref_len_dic:
        for i in range(ref_len_dic[rname]):
            consensus_pos_count_dic.setdefault(rname, {})[i + 1] = 0
    return consensus_pos_count_dic

def read_consensuspos_snv(consensuspos_snv):
    snv_count_dic = {}
    ref_consensus_pos_dic = {}
    with open(consensuspos_snv, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            ref_pos = ll[0] + "\t" + ll[1]
            consensus_pos = ll[-2], int(ll[-1])
            ref_consensus_pos_dic[ref_pos] = consensus_pos


            ref_alt_base = ll[2].upper() + ll[3].upper()
            if ref_pos in snv_count_dic and ref_alt_base in snv_count_dic[ref_pos]:
                snv_count_dic[ref_pos][ref_alt_base] += 1
            else:
                snv_count_dic.setdefault(ref_pos, {})[ref_alt_base] = 1

    return snv_count_dic, ref_consensus_pos_dic

def read_depth(depth):
    depth_dic = {}
    with open(depth, "r") as inf:
        for line in inf:
            ll = line.rstrip("\n").split("\t")
            depth_dic.setdefault(ll[0], {})[int(ll[1])] = int(ll[2])
    return depth_dic


def export_ref_pos_count_mutation_type(snv_count_dic, ref_consensus_pos_dic, consensus_pos_count_dic, depth_dic, out1):
    for ref_pos in snv_count_dic:
        mutation_count = len(snv_count_dic[ref_pos])
        consensus_chr, consensus_pos = ref_consensus_pos_dic[ref_pos]
        consensus_pos_count_dic.setdefault(consensus_chr, {})[consensus_pos] += mutation_count

    with open(out1, "w") as final1:
        for con_chr in consensus_pos_count_dic:
            for con_pos in consensus_pos_count_dic[con_chr]:
                mutation_count = consensus_pos_count_dic[con_chr][con_pos]

                if con_chr in depth_dic and con_pos in depth_dic[con_chr]:
                    depth = depth_dic[con_chr][con_pos]
                else:
                    depth = 0

                if depth != 0:
                    mutation_rate = mutation_count / depth
                else:
                    mutation_rate = 0

                if con_chr == "DMLTR5":
                    ref_ltr_name = "HMSBEAGLE"
                    ref_ltr_type = "LTR"
                else:        
                    ref_ltr_name = "_".join(re.split(r'[-_]', con_chr)[:-1])
                    ref_ltr_type = re.split(r'[-_]', con_chr)[-1]

                out_line = con_chr + "\t" + ref_ltr_name + "\t" + ref_ltr_type + "\t" + str(con_pos) + "\t" + str(mutation_count) + "\t" + str(depth) + "\t" + str(mutation_rate)
                final1.write(out_line + "\n")


if __name__ == "__main__":
    main()
