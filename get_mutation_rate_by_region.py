import sys

def main():
    snv = sys.argv[1] # AA.ds.q20.single_read.no_false_mapping.out
    total_length = sys.argv[2] # AA_AAC.genome_region.depth.qv20_qv30.length
    tag1 = sys.argv[3] # AA / AAC
    tag2 = sys.argv[4] # ds / ss
    tag3 = sys.argv[5] # q20 / q30
    dm6_bed_ls = "/rd/caiya/LTR/duplex/AAC/mismatch_py/mapq60/pileup/subreads_vaf_80/get_input_seq/bed_list"

    bed_dic = read_dm6_bed_ls(dm6_bed_ls)
    length_dic = get_total_length(total_length, tag1, tag3)
    region_count_dic = read_snv(snv, bed_dic)
    count_mutation_rate(length_dic, region_count_dic, tag1, tag2, tag3)


def read_dm6_bed_ls(dm6_bed_ls):
    bed_dic = {}
    with open(dm6_bed_ls, "r") as inf:
        for file in inf:
            file = file.rstrip("\n")
            region_type = file.split(".")[-2]
            with open(file, "r") as inf1:
                for line in inf1:
                    line = line.rstrip("\n")
                    ll = line.split("\t")
                    pos = ll[1] + "\t" + ll[2]
                    bed_dic.setdefault(ll[0], {})[pos] = region_type
    return bed_dic

def get_total_length(total_length, tag1, tag3):
    length_dic = {}
    with open(total_length, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            if ll[0] == tag1 and ll[1] == tag3:
                length_dic[ll[2]] = int(ll[3])
    return length_dic

def get_region_overlap(pos0, pos_dic):
    region_type = ""
    for pos1 in pos_dic:
        pos1_s = int(pos1.split("\t")[0])
        pos1_e = int(pos1.split("\t")[1])
        if pos0 >= pos1_s and pos0 < pos1_e:
            region_type = pos_dic[pos1]
            break
    return region_type

def read_snv(snv, bed_dic):
    region_count_dic = {}
    with open(snv, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            pos_dic = bed_dic[ll[0]]
            region_type = get_region_overlap(int(ll[1]), pos_dic)
            if region_type in region_count_dic:
                region_count_dic[region_type] += 1
            else:
                region_count_dic[region_type] = 1
    return region_count_dic

def get_mutation_rate(tag2, snv_count, region_length):
    if tag2 == "ds":
        mutation_rate = snv_count / (region_length / 2)
    elif tag2 == "ss":
        mutation_rate = snv_count / region_length
    return mutation_rate

def count_mutation_rate(length_dic, region_count_dic, tag1, tag2, tag3):
    total_snv_count = 0
    total_region_length = 0
    for region in length_dic:
        if region in region_count_dic:
            snv_count = region_count_dic[region]
        else:
            snv_count = 0
        
        region_length = length_dic[region]
        mutation_rate = get_mutation_rate(tag2, snv_count, region_length)

        total_region_length += region_length
        total_snv_count += snv_count

        out_line = tag1 + "\t" + tag2 + "\t" + tag3 + "\t" + region + "\t" + str(snv_count) + "\t" + str(region_length) + "\t" + str(mutation_rate)
        print(out_line)

    total_mutation_rate = get_mutation_rate(tag2, total_snv_count, total_region_length)
    total_out_line = tag1 + "\t" + tag2 + "\t" + tag3 + "\ttotal\t" + str(total_snv_count) + "\t" + str(total_region_length) + "\t" + str(total_mutation_rate)
    print(total_out_line)

if __name__ == "__main__":
    main()


