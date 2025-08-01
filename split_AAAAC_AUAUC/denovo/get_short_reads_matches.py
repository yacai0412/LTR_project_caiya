import sys

def main():
    insnv = sys.argv[1] # AA.indel.consensuspos.out
    short_reads_hap1_sam = sys.argv[2] # AA.indel.consensuspos.25bp.AAC.hap1.sam
    short_reads_hap2_sam = sys.argv[3] # AA.indel.consensuspos.25bp.AAC.hap2.sam

    out = ".".join(insnv.split("/")[-1].split(".")[:-1]) + ".denovo_bestmapping.out"

    hap1_mismatch_bestmapping_dic = read_sam_get_mutation_from_NMtag(short_reads_hap1_sam)
    hap2_mismatch_bestmapping_dic = read_sam_get_mutation_from_NMtag(short_reads_hap2_sam)
    filter_snv_bestmapping(insnv, hap1_mismatch_bestmapping_dic, hap2_mismatch_bestmapping_dic, out)

def read_sam_get_mutation_from_NMtag(sam):
    hap_mismatch_bestmapping_dic = {}
    with open(sam, "r") as inf:
        for line in inf:
            if not line.startswith("@"):
                ll = line.rstrip("\n").split("\t")
                snv_id = ll[0]
                cigar = ll[5]
                NMtag = ll[11].split(":")
                if int(ll[1]) == 4:
                    hap_mismatch_bestmapping_dic[snv_id] = "unmaped"
                elif "S" not in cigar and "H" not in cigar and "I" not in cigar and "D" not in cigar and int(NMtag[2]) == 0:
                    hap_mismatch_bestmapping_dic[snv_id] = "bestmapping"
                else:
                    hap_mismatch_bestmapping_dic[snv_id] = "mismatch"
    return hap_mismatch_bestmapping_dic

def filter_snv_bestmapping(insnv, hap1_mismatch_bestmapping_dic, hap2_mismatch_bestmapping_dic, out):
    with open(insnv, "r") as inf, open(out, "w") as final:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            snv_id = ll[0]
            if snv_id in hap1_mismatch_bestmapping_dic and snv_id in hap2_mismatch_bestmapping_dic:
                if  hap1_mismatch_bestmapping_dic[snv_id] == "bestmapping" or hap2_mismatch_bestmapping_dic[snv_id] == "bestmapping":
                    out_tag = "bestmapping"
                else:            
                    out_tag = "de_novo"
            else:
                print("error in " + line)


            out_line = line + "\t" + out_tag
            final.write(out_line + "\n")

if __name__ == "__main__":
    main()
