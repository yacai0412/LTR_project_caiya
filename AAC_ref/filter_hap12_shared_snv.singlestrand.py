import sys

def main():
    hap1_snv = sys.argv[1]  # AA.deepconsensus.cx3.q20.AAC.hap1.F904.overlap90.s.bam.ds.out
    hap2_snv = sys.argv[2]  # AA.deepconsensus.cx3.q20.AAC.hap2.F904.overlap90.s.bam.ds.out
    ds_ss = sys.argv[3]     # ds or ss

    out_prefix = sys.argv[4]    # AA.deepconsensus.cx3.q20.AAC.F904.overlap90.s.bam.ds

    shared_out = out_prefix + ".hap1_hap2_shared.out"    # AA.deepconsensus.cx3.q20.AAC.F904.overlap90.s.bam.ds.hap1_hap2_shared.out
    hap1_sp_out = out_prefix + ".hap1_sp.out"   # AA.deepconsensus.cx3.q20.AAC.F904.overlap90.s.bam.ds.hap1_sp.out
    hap2_sp_out = out_prefix + ".hap2_sp.out"   # AA.deepconsensus.cx3.q20.AAC.F904.overlap90.s.bam.ds.hap2_sp.out

    if ds_ss == "ds":
        hap1_read_snv_dic = read_ds_snv(hap1_snv)
        hap2_read_snv_dic = read_ds_snv(hap2_snv)
    elif ds_ss == "ss":
        hap1_read_snv_dic = read_ss_snv(hap1_snv)
        hap2_read_snv_dic = read_ss_snv(hap2_snv)

    compare_hap1_hap2_snv_dic(hap1_read_snv_dic, hap2_read_snv_dic, shared_out, hap1_sp_out, hap2_sp_out)


def get_ds_fwd_pos(plus_read, minus_read): # ll[4] ll[5]
    plus_strand = plus_read.split(":")[0].split("/")[-1]
    minus_strand = minus_read.split(":")[0].split("/")[-1]
    if plus_strand == "fwd":
        fwd_pos = ":".join(plus_read.split(":")[0:2])
    elif minus_strand == "fwd":
        fwd_pos = ":".join(minus_read.split(":")[0:2])
    else:
        print("error: " + plus_strand + "\t" + minus_read)
    return fwd_pos

def read_ds_snv(ds_snv):
    read_snv_dic = {}
    with open(ds_snv, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            fwd_pos = get_ds_fwd_pos(ll[4], ll[5])
            read_snv_dic[fwd_pos] = line
    return read_snv_dic

def read_ss_snv(ss_snv):
    read_snv_dic = {}
    with open(ss_snv, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            if ll[4] != "NA":
                ss_read_pos = ":".join(ll[4].split(":")[0:2])
            elif ll[5] != "NA":
                ss_read_pos = ":".join(ll[5].split(":")[0:2])
            read_snv_dic[ss_read_pos] = line
    return read_snv_dic

def compare_hap1_hap2_snv_dic(hap1_read_snv_dic, hap2_read_snv_dic, shared_out, hap1_sp_out, hap2_sp_out):
    with open(shared_out, "w") as shared_final, open(hap1_sp_out, "w") as hap1_sp_final, open(hap2_sp_out, "w") as hap2_sp_final:
        for read in hap1_read_snv_dic:
            hap1_line = hap1_read_snv_dic[read]
            if read in hap2_read_snv_dic:
                hap2_line = hap2_read_snv_dic[read]
                out_line = hap1_line + "\t" + hap2_line
                shared_final.write(out_line + "\n")

                hap1_read_snv_dic[read] = ""
                hap2_read_snv_dic[read] = ""
            else:
                hap1_sp_final.write(hap1_line + "\n")
                
        for read in hap2_read_snv_dic:
            hap2_line = hap2_read_snv_dic[read]
            if hap2_line != "":
                hap2_sp_final.write(hap2_line + "\n")

if __name__ == "__main__":
    main()
