import sys

def main():
    indel = sys.argv[1] # AAC.indel.hetero.vcf
    snv = sys.argv[2] # AAC.snv.hetero.vcf

    out = ".".join(snv.split("/")[-1].split(".")[:-1]) + ".inside_indel.refpos.relativepos"

    indel_pos_query_dic = read_indel_file(indel)
    read_snv_file(snv, indel_pos_query_dic, out)


def read_indel_file(indel_file):
    indel_snv_ref_pos = {}
    with open(indel_file, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            if "." not in ll[2]:
                ref_chr = ll[2].split("-")[0]
                ref_start = int(ll[2].split("-")[1])
                indel_type = ll[2].split("-")[2]
                indel_len = int(ll[2].split("-")[3])

                if indel_type == "DEL":
                    ref_end = ref_start + indel_len + 1
                elif indel_type == "INS":
                    ref_end = ref_start + 1
                else:
                    print("error in indel type" + "\t" + ll[2])

                indel_snv_ref_pos.setdefault(ref_chr, {})[(ref_start, ref_end)] = indel_type
    return indel_snv_ref_pos

def read_snv_file(snv, indel_pos_query_dic, out):
    with open(snv, "r") as inf, open(out, "w") as final:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            ref_chr = ll[2].split("-")[0]
            snv_pos = int(ll[2].split("-")[1])

            if ref_chr in indel_pos_query_dic:
                indel_dic_ref_chr = indel_pos_query_dic[ref_chr]
                min_direction, min_distace, min_indel_type = get_nearest_indel(snv_pos, indel_dic_ref_chr)
                if min_direction == "left":
                    min_distace = (-1) * min_distace
                out_line = line + "\t" + min_direction + "\t" + str(min_distace) + "\t" + str(min_indel_type)
                final.write(out_line + "\n")


def get_nearest_indel(snv_pos, indel_dic):
    min_distace = 1e100
    min_direction = ""
    min_indel_type = ""

    for indel_pos in indel_dic:
        indel_start = indel_pos[0]
        indel_end = indel_pos[1]
        indel_type = indel_dic[indel_pos]
        if snv_pos > indel_start and snv_pos < indel_end:
            min_direction = "inside"
            min_distace = min(abs(snv_pos - indel_start), abs(snv_pos - indel_end))
            min_indel_type = indel_type
            return min_direction, min_distace, min_indel_type

        elif snv_pos >= indel_start and snv_pos >= indel_end:
            direction = "right"
            distance = snv_pos - indel_end + 1
        elif snv_pos <= indel_start and snv_pos < indel_end:
            direction = "left"
            distance = indel_start - snv_pos + 1
        
        if distance <= min_distace:
            min_distace = distance
            min_direction = direction
            min_indel_type = indel_type
    return min_direction, min_distace, min_indel_type


if __name__ == "__main__":
    main()
