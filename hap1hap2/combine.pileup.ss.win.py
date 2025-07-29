import sys

def main():
    snv_file = sys.argv[1] # AA.deepconsensus.cx3.pbmm2.dm6.overlap90.s.bam.ss.qv20.out
    fwd_fwd_pileup = sys.argv[2] # AA.deepconsensus.cx3.pbmm2.dm6.overlap90.s.bam.ss.qv20.fwd_fwd.pileup
    fwd_rev_pileup = sys.argv[3] # AA.deepconsensus.cx3.pbmm2.dm6.overlap90.s.bam.ss.qv20.fwd_rev.pileup
    rev_rev_pileup = sys.argv[4] # AA.deepconsensus.cx3.pbmm2.dm6.overlap90.s.bam.ss.qv20.rev_rev.pileup
    rev_fwd_pileup = sys.argv[5] # AA.deepconsensus.cx3.pbmm2.dm6.overlap90.s.bam.ss.qv20.rev_fwd.pileup
    baseq = 93

    out = ".".join(snv_file.split("/")[-1].split(".")[:-1]) + ".plus_minus.out"
    final = open(out, "w")

    fwd_fwd_dic = read_pileup(fwd_fwd_pileup, baseq)
    fwd_rev_dic = read_pileup(fwd_rev_pileup, baseq)
    rev_rev_dic = read_pileup(rev_rev_pileup, baseq)
    rev_fwd_dic = read_pileup(rev_fwd_pileup, baseq)
    read_ss_snv(snv_file, fwd_fwd_dic, fwd_rev_dic, rev_rev_dic, rev_fwd_dic, final)
    final.close()


def get_baseq_base(base_list, baseq_list, baseq, altbase):
    out_str = ""
    for i in range(len(baseq_list)):
        bq = ord(baseq_list[i])
        base = base_list[i]
        if bq >= baseq:
            if base == "." or base == ",":
                base = altbase
            out_str = out_str + base
    return out_str

def read_pileup(pileup, baseq):
    out_dic = {}
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

def get_ss_strand(plus_reads, minus_reads):
    if plus_reads == "NA":
        support_reads_ls = minus_reads.split(",")
    else:
        support_reads_ls = plus_reads.split(",")
    return support_reads_ls

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

def get_pileup_ss(read, fwd_fwd_dic, fwd_rev_dic, rev_rev_dic, rev_fwd_dic):
    out_dic = {}
    rname = read.split(":")[0]
    rstrand = rname.split("/")[-1] # fwd/rev
    rpos = read.split(":")[1]
    refstrand = read.split(":")[2]

    if rstrand == "fwd":
        a_dic = fwd_fwd_dic
        b_dic = fwd_rev_dic
    elif rstrand == "rev":
        a_dic = rev_rev_dic
        b_dic = rev_fwd_dic
    
    if rname in a_dic and rpos in a_dic[rname] and rname in b_dic and rpos in b_dic[rname]:
        a_bases = a_dic[rname][rpos]
        a_bases1 = get_complement_base(a_bases, refstrand)
        b_bases = b_dic[rname][rpos]
        b_bases1 = get_complement_base(b_bases, refstrand)
        if a_bases1 != "" and b_bases1 != "":
            if refstrand == "+":
                out_dic["plus"] = a_bases1
                out_dic["minus"] = b_bases1
            elif refstrand == "-":
                out_dic["plus"] = b_bases1
                out_dic["minus"] = a_bases1
    return out_dic


def read_ss_snv(snv_file, fwd_fwd_dic, fwd_rev_dic, rev_rev_dic, rev_fwd_dic, final):
    with open(snv_file, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")

            support_reads_ls = get_ss_strand(ll[4], ll[5])
            plus_pileup_ls = []
            minus_pileup_ls = []
            ss_support_reads_out_ls = []
            for read in support_reads_ls:
                read_dic = get_pileup_ss(read, fwd_fwd_dic, fwd_rev_dic, rev_rev_dic, rev_fwd_dic)
                if len(read_dic) != 0:
                    plus_pileup_ls.append(read_dic["plus"])
                    minus_pileup_ls.append(read_dic["minus"])
                    ss_support_reads_out_ls.append(read)
            
            if len(plus_pileup_ls) > 0 and len(minus_pileup_ls) > 0 and len(ss_support_reads_out_ls) > 0:
                out_str = ",".join(plus_pileup_ls) + "\t" + ",".join(minus_pileup_ls) + "\t" + ",".join(ss_support_reads_out_ls)
                out_line = "\t".join(ll[0:4]) + "\t" + out_str
                final.write(out_line + "\n")

if __name__ == "__main__":
    main()
