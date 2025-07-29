import sys

def main():
    hap1_raw_snv = sys.argv[1]  # AA.deepconsensus.cx3.q20.AAC.hap1.F904.overlap90.s.bam.snv
    hap2_raw_snv = sys.argv[2]  # AA.deepconsensus.cx3.q20.AAC.hap2.F904.overlap90.s.bam.snv

    out_prefix = sys.argv[3]    # AA.deepconsensus.cx3.q20.AAC.F904.overlap90

    hap1_shared_out = out_prefix + ".shared.hap1.snv"    # AA.deepconsensus.cx3.q20.AAC.F904.overlap90.s.bam.ds.shared.hap1.snv
    hap2_shared_out = out_prefix + ".shared.hap2.snv"    # AA.deepconsensus.cx3.q20.AAC.F904.overlap90.s.bam.ds.shared.hap2.snv
    # hap1_sp_out = out_prefix + ".hap1_sp.snv"   # AA.deepconsensus.cx3.q20.AAC.F904.overlap90.s.bam.ds.hap1_sp.snv
    # hap2_sp_out = out_prefix + ".hap2_sp.snv"   # AA.deepconsensus.cx3.q20.AAC.F904.overlap90.s.bam.ds.hap2_sp.snv

    hap1_read_snv_dic = read_raw_snv(hap1_raw_snv)
    hap2_read_snv_dic = read_raw_snv(hap2_raw_snv)

    compare_hap1_hap2_snv_dic(hap1_read_snv_dic, hap2_read_snv_dic, hap1_shared_out, hap2_shared_out)


def read_raw_snv(raw_snv):
    read_snv_dic = {}
    with open(raw_snv, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            read_pos = ":".join(ll[4].split(":")[0:2])
            read_snv_dic[read_pos] = line
    return read_snv_dic

def compare_hap1_hap2_snv_dic(hap1_read_snv_dic, hap2_read_snv_dic, hap1_shared_out, hap2_shared_out):
    with open(hap1_shared_out, "w") as hap1_shared_final, open(hap2_shared_out, "w") as hap2_shared_final:
        for read in hap1_read_snv_dic:
            hap1_line = hap1_read_snv_dic[read]
            if read in hap2_read_snv_dic:
                hap2_line = hap2_read_snv_dic[read]
                # out_line = hap1_line + "\t" + hap2_line
                hap1_shared_final.write(hap1_line + "\n")
                hap2_shared_final.write(hap2_line + "\n")

                hap1_read_snv_dic[read] = ""
                hap2_read_snv_dic[read] = ""
                
if __name__ == "__main__":
    main()
