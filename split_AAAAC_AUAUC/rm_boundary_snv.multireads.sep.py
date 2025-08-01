import sys

def main():
    in_fastq = sys.argv[1:3] # /rd/caiya/LTR/duplex/multi_reads_support/supplementary/AA.deepconsensus.cx3.fwd.fastq1 /rd/caiya/LTR/duplex/multi_reads_support/supplementary/AA.deepconsensus.cx3.rev.fastq1
    sep_bed = sys.argv[3] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/ecc/AA.deepconsensus.cx3.fwd_rev.trans1_te1_ref.win.s.paf.ecc1LTR.rm_no_split.bed
    rm_length = int(sys.argv[4]) # 25?
    in_snv_list = sys.argv[5]

    sep_seq_start_0_base_dic = read_in_bed(sep_bed)
    readlength_dic = get_readlength_fastq(in_fastq)
    rm_boundary_snv(in_snv_list, sep_seq_start_0_base_dic, readlength_dic, rm_length)


def get_readlength_fastq(in_fastq):
    readlength_dic = {}
    for f in in_fastq:
        with open(f, "r") as inf:
            c = 0
            for line in inf:
                c += 1
                if c % 4 == 1:
                    line = line.rstrip("\n")
                    name = line.lstrip("@")
                    readlength_dic[name] = 0
                elif c % 4 == 2:
                    line = line.rstrip("\n")
                    readlength_dic[name] = len(line)
    return readlength_dic

def read_in_bed(inbed):
    sep_seq_start_0_base_dic = {}
    with open(inbed, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            sep_rname = ll[3]
            start = int(ll[1])
            sep_seq_start_0_base_dic[sep_rname] = start
    return sep_seq_start_0_base_dic

def rm_boundary_snv(in_snv_list, sep_seq_start_0_base_dic, readlength_dic, rm_length):
    with open(in_snv_list, "r") as snvs:
        for snv in snvs:
            snv = snv.rstrip("\n")
            out = ".".join(snv.split(".")[:-1]) + ".rm_" + str(rm_length) + "bp." + snv.split(".")[-1]
            with open(snv, "r") as inf, open(out, "w") as final:
                for line in inf:
                    line = line.rstrip("\n")
                    ll = line.split("\t")
                    plus_pileup = ll[4].split(",")
                    minus_pileup = ll[5].split(",")
                    support_reads = ll[6].split(",")

                    if len(ll) == 7: # ss DNA snv / ds DNA snv (fwd only)
                        out_plus_pileup = []
                        out_minus_pileup = []
                        out_reads = []
                        for i in range(len(support_reads)):
                            read = support_reads[i]
                            sep_rname = read.split(":")[0]
                            sep_rpos = int(read.split(":")[1])
                            whole_rname = "_".join(sep_rname.split("_")[:-1])
                            whole_rspos = sep_seq_start_0_base_dic[sep_rname] + sep_rpos

                            if whole_rname in readlength_dic:
                                readlength = readlength_dic[whole_rname]
                                relative_length = min(whole_rspos, abs(readlength - whole_rspos))
                                if relative_length > rm_length:
                                    out_plus_pileup.append(plus_pileup[i])
                                    out_minus_pileup.append(minus_pileup[i])
                                    out_reads.append(support_reads[i])
                            else:
                                print("error in\t" + snv + "\t" + line)

                        if len(out_reads) > 0:
                            out_line = "\t".join(ll[0:4]) + "\t" + ",".join(out_plus_pileup) + "\t" + ",".join(out_minus_pileup) + "\t" + ",".join(out_reads) 
                            final.write(out_line + "\n")


                    elif len(ll) == 8: # ds DNA snv
                        support_reads1 = ll[7].split(",")
                        out_plus_pileup = []
                        out_minus_pileup = []
                        out_reads = []
                        out_reads1 = []
                        for i in range(len(support_reads)):
                            read = support_reads[i]
                            sep_rname = read.split(":")[0]
                            sep_rpos = int(read.split(":")[1])
                            whole_rname = "_".join(sep_rname.split("_")[:-1])
                            whole_rspos = sep_seq_start_0_base_dic[sep_rname] + sep_rpos

                            if whole_rname in readlength_dic:
                                readlength = readlength_dic[whole_rname]
                                relative_length = min(whole_rspos, abs(readlength - whole_rspos))
                                if relative_length > rm_length:
                                    out_plus_pileup.append(plus_pileup[i])
                                    out_minus_pileup.append(minus_pileup[i])
                                    out_reads.append(support_reads[i])
                                    out_reads1.append(support_reads1[i])

                            else:
                                print("error in\t" + snv + "\t" + line)

                        if len(out_reads) > 0:
                            out_line = "\t".join(ll[0:4]) + "\t" + ",".join(out_plus_pileup) + "\t" + ",".join(out_minus_pileup) + "\t" + ",".join(out_reads) + "\t" + ",".join(out_reads1) 
                            final.write(out_line + "\n")


                    elif len(ll) == 9: # ds DNA indel (fwd only)
                        out_plus_pileup = []
                        out_minus_pileup = []
                        out_reads = []
                        out_reads1 = []
                        for i in range(len(support_reads)):
                            read = support_reads[i]
                            sep_rname = read.split(":")[0]
                            sep_rpos = int(read.split(":")[1])
                            whole_rname = "_".join(sep_rname.split("_")[:-1])
                            whole_rspos = sep_seq_start_0_base_dic[sep_rname] + sep_rpos

                            if whole_rname in readlength_dic:
                                readlength = readlength_dic[whole_rname]
                                relative_length = min(whole_rspos, abs(readlength - whole_rspos))
                                if relative_length > rm_length:
                                    final.write(line + "\n")
                            else:
                                print("error in\t" + snv + "\t" + line)


if __name__ == "__main__":
    main()

