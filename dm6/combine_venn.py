import sys

def main():
    insnv = sys.argv[1] # AA_AAC_C01_E01.ds.q20.raw.multireads.compare.out
    intag1 = sys.argv[2] # pri
    intag2 = sys.argv[3] # sup

    out = insnv + "1"
    final = open(out, "w")

    snv_dic = read_insnv(insnv)
    get_out_snv_venn(snv_dic, intag1, intag2, final)


def read_insnv(insnv):
    snv_dic = {}
    c = 0
    with open(insnv, "r") as inf:
        for line in inf:
            c += 1
            line = line.rstrip("\n")
            ll = line.split("\t")
            snvid = "snv_" + str(c)
            for tag in ll[4:]:
                snv_dic.setdefault(tag, []).append(snvid)
    return snv_dic

def get_snv_from_ls(ls, i):
    if i < len(ls):
        snv = ls[i]
    else:
        snv = ""
    return snv

def get_out_snv_venn(snv_dic, intag1, intag2, final):
    ls1 =  snv_dic[intag1]
    ls2 =  snv_dic[intag2]
    
    final.write(intag1 + "\t" + intag2 + "\n")
    mlength = max(len(ls1), len(ls2))
    for i in range(mlength):
        snv1 = get_snv_from_ls(ls1, i)
        snv2 = get_snv_from_ls(ls2, i)
        out_line = snv1 + "\t" + snv2
        final.write(out_line + "\n")


if __name__ == "__main__":
    main()
