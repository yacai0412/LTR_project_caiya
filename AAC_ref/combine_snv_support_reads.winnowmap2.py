import sys

def main():
    snv = sys.argv[1]

    out = snv.split("/")[-1] + "0"
    final = open(out, "w")

    snv_reads_dic = read_snv(snv)
    combine_snv_reads_dic(snv_reads_dic, final)
    final.close()

def read_snv(snv):
    snv_reads_dic = {}
    with open(snv, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            read = ll[4]
            strand = ll[4].split(":")[2]
            snv_pos = ll[0] + "\t" + ll[1]
            mutation_direction = ll[2].upper() + "\t" + ll[3].upper()
            snv_pos_mutation_direction = snv_pos + "\t" + mutation_direction
            if snv_pos_mutation_direction in snv_reads_dic and strand in snv_reads_dic[snv_pos_mutation_direction]:
                snv_reads_dic[snv_pos_mutation_direction][strand] = snv_reads_dic[snv_pos_mutation_direction][strand] + "," + read
            else:
                snv_reads_dic.setdefault(snv_pos_mutation_direction, {})[strand] = read
    return snv_reads_dic


def combine_snv_reads_dic(snv_reads_dic, final):
    out_dic = {}
    for snv_pos in snv_reads_dic:
        snv_chr = snv_pos.split("\t")[0]
        snv_pos0 = int(snv_pos.split("\t")[1])
        plus_reads = "NA"
        minus_reads = "NA"
        if "+" in snv_reads_dic[snv_pos]:
            plus_reads = snv_reads_dic[snv_pos]["+"]
        if "-" in snv_reads_dic[snv_pos]:
            minus_reads = snv_reads_dic[snv_pos]["-"]
        
        out_line = snv_pos + "\t" + plus_reads + "\t" + minus_reads
        if snv_chr in out_dic and snv_pos0 in out_dic[snv_chr]:
            out_dic[snv_chr][snv_pos0] = out_dic[snv_chr][snv_pos0] + "|" + out_line
        else:
            out_dic.setdefault(snv_chr, {})[snv_pos0] = out_line
    
    for snv_chr in sorted(out_dic):
        for snv_pos0 in sorted(out_dic[snv_chr]):
            out_line_ls = out_dic[snv_chr][snv_pos0].split("|")
            for out_line in out_line_ls:
                final.write(out_line + "\n")

if __name__ == "__main__":
    main()
