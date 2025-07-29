import sys

def main():
    combine_vcf = sys.argv[1] # G0.G0-F100.G73.sv.vcf

    out_vcf1 = combine_vcf.replace(".vcf", ".singleton.strict.vcf")
    out_vcf2 = combine_vcf.replace(".vcf", ".singleton.loose.vcf")

    read_vcf(combine_vcf, out_vcf1, out_vcf2)


def read_vcf(combine_vcf, out_vcf1, out_vcf2):
    with open(combine_vcf, 'r') as f, open(out_vcf1, "w") as outf1, open(out_vcf2, "w") as outf2:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith('##'):
                outf1.write(line + "\n")
                outf2.write(line + "\n")
            elif line.startswith('#CHROM'):
                header = line.strip().split('\t')
                sample_names = header[9:]
                outf1.write(line + "\tshared_specific\n")
                outf2.write(line + "\tshared_specific\n")
            else:
                line = line.rstrip("\n")
                ll = line.split("\t")
                gt_ls = ll[9:]
                
                gt_count_dic = {}
                for gt in gt_ls:
                    gt_count_dic[gt] = gt_count_dic.get(gt, 0) + 1
                
                out_tag_ls = []
                for i in range(len(sample_names)):
                    sample_name = sample_names[i]
                    gt0 = gt_ls[i]
                    
                    if gt0 != "." and gt_count_dic[gt0] == 1:
                        singleton_tag = sample_name + "_singleton"
                        out_tag_ls.append(singleton_tag)
                
                if len(out_tag_ls) == 0:
                    shared_specific_tag1 = "shared"
                    shared_specific_tag2 = "shared"
                else:
                    shared_specific_tag2 = ",".join(out_tag_ls)
                    # if "." in gt_count_dic and gt_count_dic["."] >= len(sample_names) - 1:
                    if "." in gt_count_dic:
                        shared_specific_tag1 = "unknown"
                    else:
                        shared_specific_tag1 = ",".join(out_tag_ls)

                outf1.write(line + "\t" + shared_specific_tag1 + "\n")
                outf2.write(line + "\t" + shared_specific_tag2 + "\n")


if __name__ == "__main__":
    main()
