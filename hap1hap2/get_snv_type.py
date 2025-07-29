import sys

def main():
    snv_dm6_pos = sys.argv[1] # AA.deepconsensus.cx3.q20.fq.F904.overlap90.s.bam.ds.plus_minus.subpass_5.vaf_0.8.rm_20bp.LR100.bwa.pos
    dm6_bed = "/rd/caiya/LTR/duplex/AAC_ref/dm6.total.bed"

    out = snv_dm6_pos + "1"

    dm6_bed_dic = read_dm6_bed(dm6_bed)
    read_pos_file(snv_dm6_pos, dm6_bed_dic, out)


def read_dm6_bed(dm6_bed):
    dm6_bed_dic = {}
    with open(dm6_bed, "r") as inf:
        for line in inf:
            ll = line.rstrip("\n").split("\t")
            dm6_bed_dic.setdefault(ll[0], {})[ll[1]+"\t"+ll[2]] = ll[3]
    return dm6_bed_dic


def read_pos_file(snv_dm6_pos, dm6_bed_dic, out):
    with open(snv_dm6_pos, "r") as inf, open(out, "w") as final:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            out_type = "error"
            if ll[0] in dm6_bed_dic:
                for pos in dm6_bed_dic[ll[0]]:
                    poss = int(pos.split("\t")[0])
                    pose = int(pos.split("\t")[1])
                    if int(ll[1]) >= poss and int(ll[1]) >= pose:
                        out_type = dm6_bed_dic[ll[0]][pos]
                        break

            out_line = line + "\t" + out_type
            final.write(out_line + "\n")

if __name__ == "__main__":
    main()
