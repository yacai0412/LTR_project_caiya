import sys

def main():
    inbed_path_smaplename = sys.argv[1]

    out = sys.argv[2]

    bed_files_dic, samples_ls = read_inbed_path_samplename(inbed_path_smaplename)
    total_bed_pos_dic = read_bed_files_dic(bed_files_dic)
    combine_samples_bed_pos_dic(total_bed_pos_dic, samples_ls, out)


def read_inbed_path_samplename(inbed_path_smaplename):
    bed_files_dic = {}
    samples_ls = []
    with open(inbed_path_smaplename, "r") as inf0:
        for line in inf0:
            line = line.rstrip("\n")
            ll = line.split("\t")
            sample_name = ll[0]
            bed_path = ll[1]

            bed_files_dic[sample_name] = bed_path
            samples_ls.append(sample_name)
    return bed_files_dic, samples_ls

def read_bed_files_dic(bed_files_dic):
    total_bed_pos_dic = {}
    
    for sample_name in bed_files_dic:
        bedpath = bed_files_dic[sample_name]
        with open(bedpath, "r") as inf:
            for line in inf:
                line = line.rstrip("\n")
                ll = line.split("\t")
                chrom = ll[0]
                start = ll[1]
                end = ll[2]
                te_type = ll[3]

                bedpos = chrom + "\t" + start
                total_bed_pos_dic.setdefault(bedpos, {})[sample_name] = te_type
    return total_bed_pos_dic


def combine_samples_bed_pos_dic(total_bed_pos_dic, samples_ls, out):
    with open(out, "w") as final:
        header = "#CHROM\tstart\t" + "\t".join(samples_ls) + "\tsingleton\n"
        final.write(header)
        for bedpos in total_bed_pos_dic:
            te_type_ls = []

            singleton_name = ""
            singleton_count = 0
            for sample_name in samples_ls:
                te_type = total_bed_pos_dic[bedpos].get(sample_name, "NA")
                te_type_ls.append(te_type)
                if te_type != "NA":
                    singleton_name = sample_name + ".singleton"
                    singleton_count += 1

            if singleton_count == 1:
                singleton_name = singleton_name
            else:
                singleton_name  = "shared"

            out_line = bedpos + "\t" + "\t".join(te_type_ls) + "\t" + singleton_name
            final.write(out_line + "\n")

if __name__ == "__main__":
    main()
