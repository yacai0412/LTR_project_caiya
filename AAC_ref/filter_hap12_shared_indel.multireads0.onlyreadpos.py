import sys

def main():
    hap1_raw_indel = sys.argv[1]  # AA.deepconsensus.cx3.q30.ecc1LTR1.no_filter_softclip.split.fwd.AAC.hap1.F904.s.bam.indel
    hap2_raw_indel = sys.argv[2]  # AA.deepconsensus.cx3.q30.ecc1LTR1.no_filter_softclip.split.fwd.AAC.hap2.F904.s.bam.indel

    out_prefix = sys.argv[3]    # AA.deepconsensus.cx3.q30.ecc1LTR1.no_filter_softclip.split.fwd.AAC.F904

    hap1_shared_out = out_prefix + ".shared_readpos.hap1.indel"    # $prefix.s.bam.ds.shared.hap1.indel
    hap2_shared_out = out_prefix + ".shared_readpos.hap2.indel"    # $prefix.s.bam.ds.shared.hap2.indel

    hap1_read_indel_dic = read_raw_snv_indel(hap1_raw_indel)
    hap2_read_indel_dic = read_raw_snv_indel(hap2_raw_indel)

    compare_hap1_hap2_indel_dic(hap1_read_indel_dic, hap2_read_indel_dic, hap1_shared_out, hap2_shared_out)


def read_raw_snv_indel(raw_snv_indel):
    read_snv_dic = {}
    with open(raw_snv_indel, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            ref_alt = ll[2] + "\t" + ll[3]
            read_pos = ":".join(ll[4].split(":")[0:2])
            read_snv_dic[read_pos] = line, ref_alt
    return read_snv_dic
                

def compare_hap1_hap2_indel_dic(hap1_read_indel_dic, hap2_read_indel_dic, hap1_shared_out, hap2_shared_out):
    with open(hap1_shared_out, "w") as hap1_shared_final, open(hap2_shared_out, "w") as hap2_shared_final:
        for read in hap1_read_indel_dic:
            hap1_line, hap1_ref_alt = hap1_read_indel_dic[read]
            if read in hap2_read_indel_dic:
                hap2_line, hap2_ref_alt = hap2_read_indel_dic[read]

                hap1_shared_final.write(hap1_line + "\n")
                hap2_shared_final.write(hap2_line + "\n")
                hap1_read_indel_dic[read] = ""
                hap2_read_indel_dic[read] = ""


if __name__ == "__main__":
    main()
