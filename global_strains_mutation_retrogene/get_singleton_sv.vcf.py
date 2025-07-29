import sys
import gzip

def main():
    combine_vcf = sys.argv[1] # /public3/home/a6s001736/LTR/global_strain/G0_strain_reads_LTR_ins/G0.G0_F100.G73.germline.s.vcf.gz

    # out_vcf1 = combine_vcf.replace(".vcf", ".singleton.strict.vcf")
    # out_vcf2 = combine_vcf.replace(".vcf", ".singleton.loose.vcf")
    if combine_vcf.endswith('.gz'):
        out_vcf = combine_vcf.replace(".vcf.gz", ".singleton.vcf")
    else:
        out_vcf = combine_vcf.replace(".vcf", ".singleton.vcf")

    read_vcf(combine_vcf, out_vcf)


def smart_open(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode.replace('t', ''))

def read_vcf(combine_vcf, out_vcf):
    with smart_open(combine_vcf, 'rt') as f, open(out_vcf, "w") as outf:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith('##'):
                outf.write(line + "\n")
            elif line.startswith('#CHROM'):
                header = line.strip().split('\t')
                sample_names = header[9:]
                outf.write(line + "\tshared_specific\n")
            else:
                line = line.rstrip("\n")
                ll = line.split("\t")
                gt_ls = ll[9:]
                
                gt_count_dic = {}
                for gt in gt_ls:
                    gt00 = gt.split(":")[0]
                    gt_count_dic[gt00] = gt_count_dic.get(gt00, 0) + 1
                
                out_tag_ls = []
                for i in range(len(sample_names)):
                    sample_name = sample_names[i]
                    gt0 = gt_ls[i].split(":")[0]
                    
                    if gt0 != "." and gt0 != "./." and gt_count_dic[gt0] == 1 and gt_count_dic["./."] == len(gt_ls) - 1:
                        singleton_tag = sample_name + "_singleton"
                        out_tag_ls.append(singleton_tag)
                
                if len(out_tag_ls) == 0:
                    shared_specific_tag = "shared"
                else:
                    shared_specific_tag = ",".join(out_tag_ls)

                outf.write(line + "\t" + shared_specific_tag + "\n")


if __name__ == "__main__":
    main()
