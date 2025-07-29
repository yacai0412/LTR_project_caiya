import sys

def main():
    in_af_tsv = sys.argv[1] # /public5/home/sch5655/global_strains/G0_G0_F100_G73/G0.G0-F100.G73.trEMOLO.merged.af.tsv

    out0 = ".".join(in_af_tsv.split("/")[-1].split(".")[:-2]) + ".af_over_0.tsv"
    out1 = ".".join(in_af_tsv.split("/")[-1].split(".")[:-2]) + ".af_equal_0.tsv"

    af_0_dic = read_af_tsv(in_af_tsv, out0)
    output_af_0_dic(af_0_dic, out1)

def get_af_0_outtag(TEfamily, samplename, af_0_dic, ins_pos1):
    if TEfamily == "NA":
        out_TEfamily = "NA"
    else:
        TEfamily_ls = TEfamily.split("|")
        af = float(TEfamily_ls[-1])
        if af == 0:
            af_0_dic.setdefault(samplename, {})[ins_pos1] = TEfamily
            out_TEfamily = "NA"
        else:
            out_TEfamily = TEfamily
    return af_0_dic, out_TEfamily

def read_af_tsv(in_af_tsv, out0):
    af_0_dic = {}
    with open(in_af_tsv, "r") as inf, open(out0, "w") as outf:
        for line in inf:
            if line.startswith("#"):
                outf.write(line)
            else:
                line = line.rstrip("\n")
                ll = line.split("\t")

                chrom = ll[0]
                ins_pos = ll[1]
                ins_pos1 = ll[0] + "\t" + ll[1]
                AAC_TEfamily = ll[2]
                AUC_TEfamily = ll[3]
                # G73_TEfamily = ll[4]
                # singleton_type = ll[5]

                af_0_dic, AAC_out_TEfamily = get_af_0_outtag(AAC_TEfamily, "AAC", af_0_dic, ins_pos1)
                af_0_dic, AUC_out_TEfamily = get_af_0_outtag(AUC_TEfamily, "AUC", af_0_dic, ins_pos1)
                # af_0_dic, G73_out_TEfamily = get_af_0_outtag(G73_TEfamily, "G73", af_0_dic, ins_pos1)
                out_family_dic = {"AAC": AAC_out_TEfamily, "AUC": AUC_out_TEfamily}

                if AAC_out_TEfamily == "NA" and AUC_out_TEfamily == "NA":
                    continue
                else:
                    singleton_name = ""
                    singleton_count = 0
                    for sample_name in out_family_dic:
                        te_type = out_family_dic[sample_name]
                        if te_type != "NA":
                            singleton_name = sample_name + ".singleton"
                            singleton_count += 1

                    if singleton_count == 1:
                        singleton_name = singleton_name
                    else:
                        singleton_name  = "shared"

                    outline = chrom + "\t" + ins_pos + "\t" + AAC_out_TEfamily + "\t" + AUC_out_TEfamily + "\t" + singleton_name
                    outf.write(outline + "\n")
    return af_0_dic

def output_af_0_dic(af_0_dic, out1):
    with open(out1, "w") as outf:
        for sample_name in af_0_dic:
            for ins_pos in sorted(af_0_dic[sample_name]):
                outline = ins_pos + "\t" + af_0_dic[sample_name][ins_pos] + "\t" + sample_name
                outf.write(outline + "\n")

if __name__ == "__main__":
    main()
