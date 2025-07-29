import sys

def main():
    pri_snv = sys.argv[1] # /rd/caiya/LTR/duplex/AA/mismatch/all_reads/AA.deepconsensus.cx3.pbmm2.dm6.chr.overlap90.s.bam.snv
    sup_snv = sys.argv[2] # /rd/caiya/LTR/duplex/multi_reads_support/supplementary/AA.deepconsensus.cx3.pbmm2.dm6.chr.overlap90.s.bam.snv

    snv_dic = {}
    read_snv(pri_snv, snv_dic, "pri")
    read_snv(sup_snv, snv_dic, "sup")
    print_snv_dic(snv_dic)


def read_snv(snv, snv_dic, tag):
    tmp_dic = {}
    with open(snv, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            spos = ll[0] + "\t" + ll[1]
            mutation_direction = ll[2].upper() + "\t" + ll[3].upper()
            if spos in tmp_dic and mutation_direction in tmp_dic[spos]:
                continue
            else:
                tmp_dic.setdefault(spos, {})[mutation_direction] = 1
                if spos in snv_dic and mutation_direction in snv_dic[spos]:
                    snv_dic[spos][mutation_direction] = snv_dic[spos][mutation_direction] + "\t" + tag
                else:
                    snv_dic.setdefault(spos, {})[mutation_direction] = tag

def print_snv_dic(snv_dic):
    for spos in snv_dic:
        for mutation_direction in snv_dic[spos]:
            out_line = spos + "\t" + mutation_direction + "\t" + snv_dic[spos][mutation_direction]
            print(out_line)

def write_sp_snv(snv_dic, pri_line_dic, sup_line_dic, pri_final, sup_final):
    for spos in snv_dic:
        for mutation_direction in snv_dic[spos]:
            tags_ls = snv_dic[spos][mutation_direction].split("\t")
            if len(tags_ls) == 1:
                tag = tags_ls[0]
                if tag == "pri":
                    out_line = pri_line_dic[spos][mutation_direction]
                    pri_final.write(out_line + "\n")
                elif tag == "sup":
                    out_line = sup_line_dic[spos][mutation_direction]
                    sup_final.write(out_line + "\n")


if __name__ == "__main__":
    main()

