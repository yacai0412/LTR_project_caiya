import sys

def main():
    tei_merged_tsv = sys.argv[1] # /public5/home/sch5655/global_strains/G0_G0_F100_G73/G0.G0-F100.G73.trEMOLO.merged.tsv
    dm6_bam_reads_tab_samplename_file = sys.argv[2] 

    out = ".".join(tei_merged_tsv.split("/")[-1].split(".")[:-1]) + ".af.tsv"
    out_alt_read_name = ".".join(tei_merged_tsv.split("/")[-1].split(".")[:-1]) + ".alt.readname"

    sample_reads_tab_dic, reads_tab_dic1 = read_dm6_bam_reads_tab_files(dm6_bam_reads_tab_samplename_file)
    sample_name_alt_read_name_dic = calculate_af(tei_merged_tsv, sample_reads_tab_dic, reads_tab_dic1, out)
    output_alt_readname(sample_name_alt_read_name_dic, out_alt_read_name)


def read_dm6_bam_reads_tab_files(dm6_bam_reads_tab_samplename_file):
    sample_reads_tab_dic = {}
    reads_tab_dic1 = {}
    with open(dm6_bam_reads_tab_samplename_file, "r") as inf0:
        for line in inf0:
            line = line.rstrip("\n")
            ll0 = line.split("\t")
            sample_name = ll0[0]
            dm6_bam_reads_tab = ll0[1]

            with open(dm6_bam_reads_tab, "r") as inf1:
                for line in inf1:
                    ll1 = line.rstrip("\n").split("\t")
                    ins_pos = ll1[0].split(":")[0] + ":" + ll1[0].split(":")[1]
                    read_name = ll1[0].split(":")[2]

                    reads_tab_dic1.setdefault(sample_name, {})[read_name] = ll1

                    if ins_pos in sample_reads_tab_dic and sample_name in sample_reads_tab_dic[ins_pos]:
                        sample_reads_tab_dic[ins_pos][sample_name].append(read_name)
                    else:
                        sample_reads_tab_dic.setdefault(ins_pos, {}).setdefault(sample_name, [])        
                        sample_reads_tab_dic[ins_pos][sample_name].append(read_name)
    return sample_reads_tab_dic, reads_tab_dic1


def get_if_reads_contain_TEI(read_info_ll1, G0_tefamily1):
    part_name_ls = read_info_ll1[4].split(",")
    c = 0
    for i in part_name_ls:
        if i == G0_tefamily1:
            c += 1
            return True
    if c == 0:
        return False

def calculate_af_1(sample_name, TEfamily, sample_reads_tab_dic, ins_pos1, reads_tab_dic1):
    alt_read_name_ls = []
    if TEfamily != "NA":
        alt_count = 0
        ref_count = 0
        TEfamily1 = TEfamily.split("|")[0]
        if ins_pos1 in sample_reads_tab_dic:
            for readname in sample_reads_tab_dic[ins_pos1].get(sample_name, []):
                read_info_ll1 = reads_tab_dic1[sample_name][readname]
                if get_if_reads_contain_TEI(read_info_ll1, TEfamily1):
                    alt_count += 1
                    alt_read_name_ls.append(readname)
                else:
                    ref_count += 1
        
        af = alt_count / (alt_count + ref_count) if (alt_count + ref_count) > 0 else 0
        out_tag = TEfamily + "|" + str(alt_count) + "|" + str(ref_count) + "|" + str(af)
        return out_tag, alt_read_name_ls
    else:
        return "NA", alt_read_name_ls


def calculate_af(tei_merged_tsv, sample_reads_tab_dic, reads_tab_dic1, out):
    sample_name_alt_read_name_dic = {}
    with open(tei_merged_tsv, "r") as inf, open(out, "w") as outf:
        for line in inf:
            if line.startswith("#"):
                outf.write(line)
            else:
                line = line.rstrip("\n")
                ll = line.split("\t")

                chrom = ll[0]
                ins_pos = ll[1]
                AAC_TEfamily = ll[2]
                AUC_TEfamily = ll[3]
                # G73_TEfamily = ll[4]
                singleton_type = ll[-1]

                ins_pos1 = chrom + ":" + ins_pos
                AAC_out_tag, AAC_alt_read_name_ls = calculate_af_1("AAC", AAC_TEfamily, sample_reads_tab_dic, ins_pos1, reads_tab_dic1)
                AUC_out_tag, AUC_alt_read_name_ls = calculate_af_1("AUC", AUC_TEfamily, sample_reads_tab_dic, ins_pos1, reads_tab_dic1)
                # G73_out_tag, G73_alt_read_name_ls = calculate_af_1("G73", G73_TEfamily, sample_reads_tab_dic, ins_pos1, reads_tab_dic1)
                # outline = ll[0] + "\t" + ll[1] + "\t" + G0_out_tag + "\t" + G0_F100_out_tag + "\t" + G73_out_tag + "\t" + singleton_type
                outline = ll[0] + "\t" + ll[1] + "\t" + AAC_out_tag + "\t" + AUC_out_tag + "\t" + singleton_type
                outf.write(outline + "\n")

                sample_name_alt_read_name_dic.setdefault(ins_pos1, {})["AAC"] = AAC_alt_read_name_ls
                sample_name_alt_read_name_dic.setdefault(ins_pos1, {})["AUC"] = AUC_alt_read_name_ls
                # sample_name_alt_read_name_dic.setdefault(ins_pos1, {})["G73"] = G73_alt_read_name_ls

    return sample_name_alt_read_name_dic

def output_alt_readname(sample_name_alt_read_name_dic, out_alt_read_name):
    with open(out_alt_read_name, "w") as outf:
        for ins_pos in sample_name_alt_read_name_dic:
            for samplename in sample_name_alt_read_name_dic[ins_pos]:
                alt_read_name_ls = sample_name_alt_read_name_dic[ins_pos][samplename]
                if len(alt_read_name_ls) > 0:
                    outline = ins_pos + "\t" + samplename + "\t" + str(len(alt_read_name_ls)) + "\t" + ",".join(alt_read_name_ls)
                    outf.write(outline + "\n")


if __name__ == "__main__":
    main()
