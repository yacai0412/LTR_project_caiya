import sys

def main():
    hap1_tab = sys.argv[1] # AA.ori.ccs.fq.AAC.hap1.win.s.uniq.paf.asm.combine.tab
    hap2_tab = sys.argv[2] # AA.ori.ccs.fq.AAC.hap2.win.s.uniq.paf.asm.combine.tab

    hap1_out = ".".join(hap1_tab.split("/")[-1].split(".")[:-1]) + ".shared.chimeric.tab"
    hap2_out = ".".join(hap2_tab.split("/")[-1].split(".")[:-1]) + ".shared.chimeric.tab"

    hap1_non_chimeric_out = ".".join(hap1_tab.split("/")[-1].split(".")[:-1]) + ".non_chimeric.tab"
    hap2_non_chimeric_out = ".".join(hap2_tab.split("/")[-1].split(".")[:-1]) + ".non_chimeric.tab"


    hap1_tab_dic = read_intab(hap1_tab)
    hap2_tab_dic = read_intab(hap2_tab)

    get_hap1_hap2_shared_chimeric_reads(hap1_tab_dic, hap2_tab_dic, hap1_non_chimeric_out, hap2_non_chimeric_out, hap1_out, hap2_out)


def read_intab(intab):
    tab_chimeric_dic = {}
    with open(intab, "r") as inf:
        for line in inf:
            chimeric_tag = ""
            line = line.rstrip("\n")
            ll = line.split("\t")

            rname = ll[0] + "\t" + ll[1]
            query_pos_ls = ll[2].split(",")

            if len(query_pos_ls) == 1:
                chimeric_tag = False
            else:
                chimeric_tag = True

            tab_chimeric_dic[rname] = chimeric_tag, line
    return tab_chimeric_dic

def get_hap1_hap2_shared_chimeric_reads(hap1_tab_dic, hap2_tab_dic, hap1_non_chimeric_out, hap2_non_chimeric_out, hap1_out, hap2_out):
    with open(hap1_out, "w") as hap1_outf, open(hap2_out, "w") as hap2_outf, open(hap1_non_chimeric_out, "w") as hap1_nonchimeric_outf, open(hap2_non_chimeric_out, "w") as hap2_nonchimeric_outf:
        for rname in hap1_tab_dic:
            if rname in hap2_tab_dic and hap1_tab_dic[rname][0] and hap2_tab_dic[rname][0]:
                hap1_outf.write(hap1_tab_dic[rname][1] + "\n")
                hap2_outf.write(hap2_tab_dic[rname][1] + "\n")
            
            elif not hap1_tab_dic[rname][0]:
                hap1_nonchimeric_outf.write(hap1_tab_dic[rname][1] + "\n")

        for rname in hap2_tab_dic:
            if not hap2_tab_dic[rname][0]:
                hap2_nonchimeric_outf.write(hap2_tab_dic[rname][1] + "\n")


if __name__ == "__main__":
    main()

