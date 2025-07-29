import sys

def main():
    snv_pileup = sys.argv[1] # AA.deepconsensus.cx3.pbmm2.dm6.overlap90.s.bam.ds.qv20.plus_minus.out
    ds_ss_tag = sys.argv[2] # ds / ss
    support_subreads = int(sys.argv[3]) # 5
    subreads_vaf = float(sys.argv[4]) # 0.8

    out = ".".join(snv_pileup.split("/")[-1].split(".")[:-1]) + ".subpass_" + str(support_subreads) + ".vaf_" + str(subreads_vaf) + ".out"
    final = open(out, "w")

    if ds_ss_tag == "ds":
        read_ds_snv(snv_pileup, support_subreads, subreads_vaf, final)
    elif ds_ss_tag == "ss":
        read_ss_snv(snv_pileup, support_subreads, subreads_vaf, final)
    final.close()


def if_surbeads_pass(subreads_str, target_base, support_subreads, subreads_vaf):
    plus_subreads_ls_i = list(subreads_str)
    c_alt = 0
    c_total = 0
    for j in plus_subreads_ls_i:
        c_total += 1
        if j.upper() == target_base:
            c_alt += 1
    if c_alt >= support_subreads and c_alt/c_total >= subreads_vaf:
        return True
    else:
        return False

def read_ds_snv(snv_pileup, support_subreads, subreads_vaf, final):
    with open(snv_pileup, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            alt_base = ll[3].upper()

            plus_subreads_ls = ll[4].split(",")
            minus_subreads_ls = ll[5].split(",")
            plus_reads_ls = ll[6].split(",")
            minus_reads_ls = ll[7].split(",")

            out_plus_subreads_ls = []
            out_minus_subreads_ls = []
            out_plus_reads_ls = []
            out_minus_reads_ls =[]
            if len(plus_subreads_ls) == len(minus_subreads_ls) and len(plus_subreads_ls) == len(plus_reads_ls) and len(plus_subreads_ls) == len(minus_reads_ls):
                for i in range(len(plus_subreads_ls)):
                    if if_surbeads_pass(plus_subreads_ls[i], alt_base, support_subreads, subreads_vaf) and if_surbeads_pass(minus_subreads_ls[i], alt_base, support_subreads, subreads_vaf):
                        out_plus_subreads_ls.append(plus_subreads_ls[i])
                        out_minus_subreads_ls.append(minus_subreads_ls[i])
                        out_plus_reads_ls.append(plus_reads_ls[i])
                        out_minus_reads_ls.append(minus_reads_ls[i])
            else:
                print("error in " + line)

            if len(out_plus_subreads_ls) != 0 and len(out_minus_subreads_ls) != 0 and len(out_plus_reads_ls) != 0 and len(out_minus_reads_ls) != 0:
                out_str = ",".join(out_plus_subreads_ls) + "\t" + ",".join(out_minus_subreads_ls) + "\t" + ",".join(out_plus_reads_ls) + "\t" + ",".join(out_minus_reads_ls)
                out_line = "\t".join(ll[0:4]) + "\t" + out_str
                final.write(out_line + "\n")

def read_ss_snv(snv_pileup, support_subreads, subreads_vaf, final):
    with open(snv_pileup, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            ref_base = ll[2].upper()
            alt_base = ll[3].upper()

            plus_subreads_ls = ll[4].split(",")
            minus_subreads_ls = ll[5].split(",")
            support_reads_ls = ll[6].split(",")
            ss_ref_strand = support_reads_ls[0].split(":")[-1]

            out_plus_subreads_ls = []
            out_minus_subreads_ls = []
            out_support_reads_ls = []
            if len(plus_subreads_ls) == len(minus_subreads_ls) and len(plus_subreads_ls) == len(support_reads_ls):
                for i in range(len(plus_subreads_ls)):
                    if (ss_ref_strand == "+" and if_surbeads_pass(plus_subreads_ls[i], alt_base, support_subreads, subreads_vaf) and if_surbeads_pass(minus_subreads_ls[i], ref_base, support_subreads, subreads_vaf)) or (ss_ref_strand == "-" and if_surbeads_pass(plus_subreads_ls[i], ref_base, support_subreads, subreads_vaf) and if_surbeads_pass(minus_subreads_ls[i], alt_base, support_subreads, subreads_vaf)):
                        out_plus_subreads_ls.append(plus_subreads_ls[i])
                        out_minus_subreads_ls.append(minus_subreads_ls[i])
                        out_support_reads_ls.append(support_reads_ls[i])
            else:
                print("error in " + line)

            if len(out_plus_subreads_ls) != 0 and len(out_minus_subreads_ls) != 0 and len(out_support_reads_ls) != 0:
                out_str = ",".join(out_plus_subreads_ls) + "\t" + ",".join(out_minus_subreads_ls) + "\t" + ",".join(out_support_reads_ls)
                out_line = "\t".join(ll[0:4]) + "\t" + out_str
                final.write(out_line + "\n")

if __name__ == "__main__":
    main()
