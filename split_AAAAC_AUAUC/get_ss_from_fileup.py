import sys
import gzip


def main():
    inpileup = sys.argv[1] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/ecc/consensus/AA.deepconsensus.cx3.q30.ecc1LTR.chimeric.split.melTE.splice.k20.F904.s.bam.pileup.gz
    
    out = ".".join(inpileup.split("/")[-1].split(".")[:-1]) + ".ss.snv"

    
    read_inpileup(inpileup, out)


def reverse_complement(seq):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement_seq = [complement_dict[nuc] for nuc in seq]
    rev_comp_seq = ''.join(complement_seq[::-1])
    return rev_comp_seq

def read_inpileup(inpileup, out):
    with gzip.open(inpileup, 'rt') as inf, open(out, "w") as final:
        for line in inf:
            tmp_dic = {}

            ll = line.rstrip("\n").split("\t")
            chrom = ll[0]
            ref_pos = ll[1]
            ref_base = ll[2].upper()
            if ref_base == "A" or ref_base == "T" or ref_base == "G" or ref_base == "C":
                # out_line = chrom + "\t" + ref_pos + "\t" + ref_base
                depth = int(ll[3])
                pileup_ls = list(ll[4])
                readname_ls = ll[6].split(",")
                if depth > 0:
                    for i in range(depth):
                        rbase = pileup_ls[i]
                        rname = readname_ls[i]
                        zmw = rname.split("/")[1]
                        rstrand = rname.split("/")[3].split("_")[0]
                        if rbase == ".":
                            rbase1 = ref_base
                            tmp_dic.setdefault(zmw, {})[rstrand] = rbase1, rname, "+"
                        elif rbase == ",":
                            rbase1 = reverse_complement(ref_base)
                            tmp_dic.setdefault(zmw, {})[rstrand] = rbase1, rname, "-"
                        elif rbase.isalpha():
                            if rbase.isupper():
                                rbase1 = rbase
                                tmp_dic.setdefault(zmw, {})[rstrand] = rbase1, rname, "+"
                            else:
                                rbase1 = reverse_complement(rbase.upper())
                                tmp_dic.setdefault(zmw, {})[rstrand] = rbase1, rname, "-"

                    for zmw in tmp_dic:
                        if len(tmp_dic[zmw]) == 2 and "fwd" in tmp_dic[zmw] and "rev" in tmp_dic[zmw]:
                            fwd_base, fwd_name, fwd_ref_strand = tmp_dic[zmw]["fwd"]
                            rev_base, rev_name, rev_ref_strand = tmp_dic[zmw]["rev"]
                            if fwd_base != reverse_complement(rev_base):
                                if fwd_ref_strand == "+" and rev_ref_strand == "-":
                                    if fwd_base == ref_base or rev_base == reverse_complement(ref_base):
                                        out_ss_tag = "0"
                                    else:
                                        out_ss_tag = "1"
                                    out_line = chrom + "\t" + ref_pos + "\t" + ref_base + "\t" + fwd_name + ":" + fwd_ref_strand + ":" + fwd_base + "\t" + rev_name + ":" + rev_ref_strand + ":" + rev_base + "\t" + out_ss_tag
                                    final.write(out_line + "\n")
                                elif rev_ref_strand == "+":
                                    if fwd_base == reverse_complement(ref_base) or rev_base == ref_base:
                                        out_ss_tag = "0"
                                    else:
                                        out_ss_tag = "1"
                                    out_line = chrom + "\t" + ref_pos + "\t" + ref_base + "\t" + rev_name + ":" + rev_ref_strand + ":" + rev_base + "\t" + fwd_name + ":" + fwd_ref_strand + ":" + fwd_base + "\t" + out_ss_tag
                                    final.write(out_line + "\n")
                                else:
                                    print("error in " + out_line + " " + fwd_name)


if __name__ == "__main__":
    main()
