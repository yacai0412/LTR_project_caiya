import sys
import pysam

def main():
    AAC_AUC_vcf_tab = sys.argv[1] # AAC_AUC.ccs.deepconsensus.asm.split.fq.melTE.mm2sr.s.paf.combine.tab1.TEI.vcf
    hap1_bam = sys.argv[2]
    hap2_bam = sys.argv[3]
    intag = sys.argv[4]  # AAC / AUC

    shared_out = ".".join(AAC_AUC_vcf_tab.split("/")[-1].split(".")[:-1]) + ".shared.vcf"
    specific_out = ".".join(AAC_AUC_vcf_tab.split("/")[-1].split(".")[:-1]) + "." + intag + ".specific.vcf"

    tei_rnames_dic, specific_reads_dic = read_AAC_AUC_vcf_tab(AAC_AUC_vcf_tab, intag)
    hap1_bam_soft_clip_dic = read_hap_bam(hap1_bam, specific_reads_dic)
    hap2_bam_soft_clip_dic = read_hap_bam(hap2_bam, specific_reads_dic)
    get_specific_insertion_segregate(hap1_bam_soft_clip_dic, hap2_bam_soft_clip_dic, tei_rnames_dic, intag, shared_out, specific_out)



def read_AAC_AUC_vcf_tab(AAC_AUC_vcf_tab, intag):
    tei_rnames_dic = {}
    specific_reads_dic = {}
    with open(AAC_AUC_vcf_tab, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")

            tei_pos_type = ll[0] + "\t" + ll[1] + "\t" + ll[2]
            tei_rnames_dic[tei_pos_type] = ll

            if intag == "AAC" and ll[4] == "NA":
                for r in ll[3].split(","):
                    specific_reads_dic[r] = intag

            if intag == "AUC" and ll[3] == "NA":
                for r in ll[4].split(","):
                    specific_reads_dic[r] = intag

    return tei_rnames_dic, specific_reads_dic

def read_hap_bam(hap_bam, specific_reads_dic):
    hap_bam_soft_clip_dic = {}
    with pysam.AlignmentFile(hap_bam, "rb") as samfile:
        for read in samfile:
            if read.query_name in specific_reads_dic and not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                softclip_tag = "not_clipped"
                for cigar_type, cigar_length in read.cigartuples:
                    if (cigar_type == 4 or cigar_type == 5) and cigar_length >= 10 :  # Soft clip / Hard clip (clipped sequence present in SEQ)
                        softclip_tag = "clipped"
                        break
                hap_bam_soft_clip_dic[read.query_name] = softclip_tag
    return hap_bam_soft_clip_dic

def get_specific_insertion_segregate(hap1_bam_soft_clip_dic, hap2_bam_soft_clip_dic, tei_rnames_dic, intag, shared_out, specific_out):
    with open(shared_out, "w") as shared_outf, open(specific_out, "w") as specific_outf:
        for tei_pos_type in tei_rnames_dic:
            ll = tei_rnames_dic[tei_pos_type]
            if ll[3] != "NA" and ll[4] != "NA":
                outline = "\t".join(ll) + "\tshared\n"
                shared_outf.write(outline)
            elif intag == "AAC" and ll[4] == "NA":
                c = 0
                for r in ll[3].split(","):
                    if r in hap1_bam_soft_clip_dic or r in hap2_bam_soft_clip_dic:
                        if hap1_bam_soft_clip_dic[r] == "not_clipped" or hap2_bam_soft_clip_dic[r] == "not_clipped":
                            outline = "\t".join(ll) + "\tfixed\n"
                            specific_outf.write(outline)
                            c += 1    
                            break
                        else:
                            outline = "\t".join(ll) + "\tsegregate\n"
                            specific_outf.write(outline)
                            c += 1
                            break

            elif intag == "AUC" and ll[3] == "NA":
                c = 0
                for r in ll[4].split(","):
                    if r in hap1_bam_soft_clip_dic or r in hap2_bam_soft_clip_dic:
                        if hap1_bam_soft_clip_dic[r] == "not_clipped" or hap2_bam_soft_clip_dic[r] == "not_clipped":
                            outline = "\t".join(ll) + "\tfixed\n"
                            specific_outf.write(outline)
                            c += 1    
                            break
                        else:
                            outline = "\t".join(ll) + "\tsegregate\n"
                            specific_outf.write(outline)
                            c += 1
                            break

                if c == 0:
                    outline = "\t".join(ll) + "\tsegregate\n"
                    specific_outf.write(outline)

if __name__ == "__main__":
    main()
