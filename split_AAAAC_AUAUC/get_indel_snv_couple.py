import sys

def main():
    indel = sys.argv[1] # AA.deepconsensus.cx3.q30.ecc1LTR1.no_filter_softclip.split.fwd.AAC.F904.shared.hap1.indel.fwd_plus_minus.subpass_8.vaf_0.8.rm_25bp.1bp_indel.consensuspos.out
    snv = sys.argv[2] # AA.deepconsensus.cx3.q30.ecc1LTR1.no_filter_softclip.split.fwd.AAC.F904.shared.hap1.fwd_plus_minus.subpass_5.vaf_0.8.rm_25bp.consensuspos.out

    out = ".".join(snv.split("/")[-1].split(".")[:-1]) + ".inside_indel.relativepos"

    indel_pos_query_dic = read_indel_file(indel)
    read_snv_file(snv, indel_pos_query_dic, out)


def read_indel_file(indel_file):
    indel_snv_query_pos = {}
    with open(indel_file, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            rname = ll[4].split(":")[0]
            indel_len = int(ll[5])
            indel_type = ll[6]

            query_start = int(ll[4].split(":")[1])
            if indel_type == "DEL":
                query_end = query_start + 1
            elif indel_type == "INS":
                query_end = query_start + indel_len + 1

            indel_snv_query_pos.setdefault(rname, {})[(query_start, query_end)] = indel_type
    return indel_snv_query_pos

def read_snv_file(snv, indel_pos_query_dic, out):
    with open(snv, "r") as inf, open(out, "w") as final:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            rname = ll[4].split(":")[0]
            snv_pos = int(ll[4].split(":")[1])
            if rname in indel_pos_query_dic:
                indel_dic_rname = indel_pos_query_dic[rname]
                min_direction, min_distace, min_indel_type = get_nearest_indel(snv_pos, indel_dic_rname)
                if min_direction == "left":
                    min_distace = (-1) * min_distace
                out_line = line + "\t" + min_direction + "\t" + str(min_distace) + "\t" + str(min_indel_type)
            # new added in 2024.11.25
            else:
                out_line = line + "\tno_indel\tNA\tNA"
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
