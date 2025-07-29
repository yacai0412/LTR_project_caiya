import sys
import pysam

def main():
    inccs_readname = sys.argv[1] # /rd/caiya/LTR/duplex/AA/mismatch/pileup/TE_pos/AA.readname
    ccs_subread_bam = sys.argv[2]
    deepconsensus_fwd_fwd_subreads_bam = sys.argv[3]

    out = inccs_readname + ".ccs_fwd.strand.subreads_bam.out"
    final = open(out, "w")

    zmw_dic = read_inccs_readname(inccs_readname)
    fwd_fwd_dic = read_deepconsensus_bam(deepconsensus_fwd_fwd_subreads_bam, zmw_dic)
    read_ccs_bam(ccs_subread_bam, zmw_dic, fwd_fwd_dic, final)
    final.close()


def read_inccs_readname(inccs_readname):
    zmw_dic = {}
    with open(inccs_readname, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            rzmw = line.split("/")[1]
            zmw_dic[rzmw] = 1
    return zmw_dic

def get_strand(read):
    if read.is_forward:
        strand = "+"
    else:
        strand = "-"
    return strand

def read_deepconsensus_bam(bamfile, zmw_dic):
    fwd_fwd_dic = {}
    inbam = pysam.AlignmentFile(bamfile, "rb")
    for read in inbam:
        qname = read.query_name
        qzmw = qname.split("/")[1]
        qsubreadspass = qname.split("/")[2]
        strand = get_strand(read)

        if qzmw in zmw_dic and qzmw not in fwd_fwd_dic:
            fwd_fwd_dic.setdefault(qzmw, {})[qsubreadspass] = strand

    return fwd_fwd_dic

def read_ccs_bam(bamfile, zmw_dic, fwd_fwd_dic, final):
    out_dic = {}
    inbam = pysam.AlignmentFile(bamfile, "rb")
    for read in inbam:
        qname = read.query_name
        qzmw = qname.split("/")[1]
        qsubreadspass = qname.split("/")[2]
        strand = get_strand(read)

        if qzmw in zmw_dic and qzmw in fwd_fwd_dic and qsubreadspass in fwd_fwd_dic[qzmw]:
            fwd_strand = fwd_fwd_dic[qzmw][qsubreadspass]
            if strand == fwd_strand:
                out_dic[qzmw] = "fwd"
            else:
                out_dic[qzmw] = "rev"
    
    for qzmw in zmw_dic:
        if qzmw in out_dic:
            out_strand = out_dic[qzmw]
        else:
            out_strand = "error"
        final.write(qzmw + "\t" + out_strand + "\n")

if __name__ == "__main__":
    main()
