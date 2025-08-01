import sys

def main():
    snv_file = sys.argv[1] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/ecc/AA.deepconsensus.cx3.q30.ecc1LTR.chimeric.split.AAC.F904.shared.hap1.ds.out
    sep_readpos0 = sys.argv[2] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/ecc/AA.deepconsensus.cx3.q30.ecc1LTR.chimeric.split.AAC.F904.shared.readpos0
    pileup_files = sys.argv[3:] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/ecc/AA.deepconsensus.cx3.q30.ecc1LTR.chimeric.split.AAC.F904.shared.fwd_fwd.pileup /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/ecc/AA.deepconsensus.cx3.q30.ecc1LTR.chimeric.split.AAC.F904.shared.rev_rev.pileup
    baseq = 93

    out = ".".join(snv_file.split("/")[-1].split(".")[:-1]) + ".plus_minus.out"

    sep_name_dic = read_readpos0(sep_readpos0)
    pileup_dic = read_pileup(pileup_files, baseq)
    read_sep_ds_snv(snv_file, pileup_dic, sep_name_dic, out)


def read_readpos0(sep_readpos0):
    sep_name_dic = {}
    with open(sep_readpos0, "r") as inf:
        for line in inf:
            ll = line.rstrip("\n").split("\t")
            sepname = ll[0]
            sep_pos = ll[1]
            whole_name = ll[2]
            whole_rpos = ll[3]
            sep_name_dic.setdefault(sepname, {})[sep_pos] = whole_name, whole_rpos
    return sep_name_dic

def get_baseq_base(base_list, baseq_list, baseq, altbase):
    out_str = ""
    for i in range(len(baseq_list)):
        bq = ord(baseq_list[i])
        base = base_list[i]
        if bq >= baseq:
            if base == ".":
                base = altbase
            out_str = out_str + base
    return out_str

def read_pileup(pileup_files, baseq):
    out_dic = {}
    for pileup in pileup_files:
        with open(pileup, "r") as inf:
            for line in inf:
                line = line.rstrip("\n")
                ll = line.split("\t")
                rname = ll[0]
                rpos = ll[1]
                altbase = ll[2].upper()
                base_list = list(ll[4].upper())
                baseq_list = list(ll[5])
                bases = get_baseq_base(base_list, baseq_list, baseq, altbase)
                out_dic.setdefault(rname, {})[rpos] = bases
    return out_dic

def get_complement_base(inbases, refstrand):
    outbases = ""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    if refstrand == "+":
        outbases = inbases
    elif refstrand == "-":
        inbases_list = list(inbases)
        for b in inbases_list:
            b = b.upper()
            if b in complement:
                outbases = outbases + complement[b]
            else:
                outbases = outbases + b
    return outbases

def get_pileup(read, sep_name_dic, pileup_dic):
    sep_rname = read.split(":")[0]
    sep_rpos = read.split(":")[1]
    refstrand = read.split(":")[2]
    if sep_rname in sep_name_dic and sep_rpos in sep_name_dic[sep_rname]:
        whole_rname, whole_rpos = sep_name_dic[sep_rname][sep_rpos]

        bases1 = ""
        if whole_rname in pileup_dic and whole_rpos in pileup_dic[whole_rname]:
            bases = pileup_dic[whole_rname][whole_rpos]
            bases1 = get_complement_base(bases, refstrand)
    return bases1

def get_zmw_bases_dic(read_ls, sep_name_dic, pileup_dic):
    zmw_dic = {}
    for read in read_ls:
        zmw = read.split("/")[1]
        bases = get_pileup(read, sep_name_dic, pileup_dic)
        if bases != "":
            zmw_dic[zmw] = bases + "\t" + read
    return zmw_dic

def get_out_bases_reads_str(plus_zmw_bases_dic, minus_zmw_bases_dic):
    out_plus_bases = []
    out_plus_reads = []
    out_minus_bases = []
    out_minus_reads = []
    for zmw in plus_zmw_bases_dic:
        if zmw in minus_zmw_bases_dic:
            plus_bases = plus_zmw_bases_dic[zmw].split("\t")[0]
            plus_read = plus_zmw_bases_dic[zmw].split("\t")[1]
            minus_bases = minus_zmw_bases_dic[zmw].split("\t")[0]
            minus_read = minus_zmw_bases_dic[zmw].split("\t")[1]
            out_plus_bases.append(plus_bases)
            out_plus_reads.append(plus_read)
            out_minus_bases.append(minus_bases)
            out_minus_reads.append(minus_read)
    if len(out_plus_bases) > 0 and len(out_plus_reads) > 0 and len(out_minus_bases) > 0 and len(out_minus_reads) > 0:
        out_plus_bases_str = ",".join(out_plus_bases)
        out_plus_reads_str = ",".join(out_plus_reads)
        out_minus_bases_str = ",".join(out_minus_bases)
        out_minus_reads_str = ",".join(out_minus_reads)
        out_str = out_plus_bases_str + "\t" + out_minus_bases_str + "\t" + out_plus_reads_str + "\t" + out_minus_reads_str
    else:
        out_str = ""
    return out_str


def read_sep_ds_snv(snv_file, pileup_dic, sep_name_dic, out):
    with open(snv_file, "r") as inf, open(out, "w") as final:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            
            plus_ls = ll[4].split(",")
            minus_ls = ll[5].split(",")
            plus_zmw_bases_dic = get_zmw_bases_dic(plus_ls, sep_name_dic, pileup_dic)
            minus_zmw_bases_dic = get_zmw_bases_dic(minus_ls, sep_name_dic, pileup_dic)

            out_str = get_out_bases_reads_str(plus_zmw_bases_dic, minus_zmw_bases_dic)
            if out_str != "":
                out_line = "\t".join(ll[0:4]) + "\t" + out_str
                final.write(out_line + "\n")


if __name__ == "__main__":
    main()
