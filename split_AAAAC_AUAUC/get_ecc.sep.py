import sys
import re

def main():
    intab = sys.argv[1] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/AA.deepconsensus.cx3.fwd_rev.trans1_te1_ref.win.s.paf.tab

    out = ".".join(intab.split("/")[-1].split(".")[:-1]) + ".ecc1LTR.tab"
    out1 = ".".join(intab.split("/")[-1].split(".")[:-1]) + ".ecc2LTR.tab"

    read_intab_get_ecc_1LTR(intab, out, out1)



def read_intab_get_ecc_1LTR(intab, out, out1):
    with open(intab, "r") as inf, open(out, "w") as final, open(out1, "w") as final1:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            target_ls = ll[4].split(",")
            strand_ls = ll[7].split(",")
            type_ls = ll[8].split(",")
            if len(target_ls) == 3 and type_ls == ["TE", "TE", "TE"] and (strand_ls == ["+", "+", "+"] or strand_ls == ["-", "-", "-"]):
                if target_ls[0] == "HMSBEAGLE_I" and target_ls[2] == "HMSBEAGLE_I" and target_ls[1] == "DMLTR5":
                    final.write(line + "\n")
                else:
                    target_0 = re.split(r'[-_]', target_ls[0])
                    target_1 = re.split(r'[-_]', target_ls[1])
                    target_2 = re.split(r'[-_]', target_ls[2])
                    if target_0[-1] == "I" and target_2[-1] == "I" and target_1[-1] == "LTR" and "_".join(target_0[:-1]) == "_".join(target_1[:-1]) and "_".join(target_0[:-1]) == "_".join(target_2[:-1]):
                        final.write(line + "\n")

            elif len(target_ls) == 4 and type_ls == ["TE", "TE", "TE", "TE"] and (strand_ls == ["+", "+", "+", "+"] or strand_ls == ["-", "-", "-", "-"]):
                if target_ls[0] == "HMSBEAGLE_I" and target_ls[3] == "HMSBEAGLE_I" and target_ls[1] == "DMLTR5" and target_ls[2] == "DMLTR5":
                    final1.write(line + "\n")
                else:
                    target_0 = re.split(r'[-_]', target_ls[0])
                    target_1 = re.split(r'[-_]', target_ls[1])
                    target_2 = re.split(r'[-_]', target_ls[2])
                    target_3 = re.split(r'[-_]', target_ls[3])
                    if target_0[-1] == "I" and target_3[-1] == "I" and (target_1[-1] == "LTR" and "_".join(target_0[:-1]) == "_".join(target_1[:-1])) and (target_2[-1] == "LTR" and "_".join(target_0[:-1]) == "_".join(target_2[:-1])) and "_".join(target_0[:-1]) == "_".join(target_3[:-1]):
                        final1.write(line + "\n")


if __name__ == "__main__":
    main()
