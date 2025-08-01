import sys
import re

def main():
    intab = sys.argv[1] # /rd2/caiya/LTR/split_AAAAC_AUAUC/gene_te_ref/AA.deepconsensus.cx3.fwd_rev.trans1_te1_ref.win.s.paf.tab

    out = ".".join(intab.split("/")[-1].split(".")[:-1]) + ".twoLTR.tab"
    out_linear = ".".join(intab.split("/")[-1].split(".")[:-1]) + ".linear2LTR.tab"
    # out_ecc = ".".join(intab.split("/")[-1].split(".")[:-1]) + ".ecc2LTR.tab"

    read_intab_get_ecc_1LTR(intab, out, out_linear)


def read_intab_get_ecc_1LTR(intab, out, out_linear):
    with open(intab, "r") as inf, open(out, "w") as final, open(out_linear, "w") as final_linear:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            target_ls = ll[4].split(",")
            if target_ls[0] == "DMLTR5" and target_ls[-1] == "DMLTR5" and "HMSBEAGLE_I" in target_ls:
                final.write(line + "\n")
                if len(set(target_ls[1:-1])) == 1:
                    final_linear.write(line + "\n")

            else:
                target_i_ls = []
                target_i_dic = {}
                for i in range(len(target_ls)):
                    target_i = re.split(r'[-_]', target_ls[i])
                    target_i_name = "_".join(target_i[:-1])
                    target_i_part = target_i[-1]
                    target_i_ls.append((target_i_name, target_i_part))
                    target_i_dic.setdefault(target_i_part, {})[target_i_name] = 1
                
                if target_i_ls[0][1] == "LTR" and target_i_ls[-1][1] == "LTR" and target_i_ls[0][0] == target_i_ls[-1][0] and "I" in target_i_dic and target_i_ls[0][0] in target_i_dic["I"]:
                    final.write(line + "\n")
                    if len(set(target_i_ls[1:-1])) == 1:
                        final_linear.write(line + "\n")



if __name__ == "__main__":
    main()
