import sys

def main():
    read_specific_snv = sys.argv[1] # SV.0328.CE2.INS.A1.alt.consensus.asm.s.bam.read_specific.snv 
    inside_indel_pos = sys.argv[2] # SV.0328.CE2.INS.A1.alt.consensus.s.readpos.out

    out = read_specific_snv + ".inside_indel.relativepos"

    indel_pos_dic = read_inside_indel(inside_indel_pos)
    read_snv(read_specific_snv, indel_pos_dic, out)


def get_indel_start_end_pos(relative_pos_column):
    relative_pos = int(relative_pos_column.split(":")[0])
    length = int(relative_pos_column.split(":")[2].split("bp_")[0])
    indel_type = relative_pos_column.split(":")[2].split("bp_")[1]

    if indel_type == "deletion":
        indel_start = relative_pos
        indel_end = relative_pos + 1
    elif indel_type == "insertion":
        indel_start = relative_pos
        indel_end = relative_pos + length - 1 

    return str(indel_start) + "\t" + str(indel_end)

def read_inside_indel(inside_indel_pos):
    indel_pos_dic = {}
    with open(inside_indel_pos, "r") as inf:
        for line in inf:
            ll = line.rstrip("\n").split("\t")
            rname = ll[0]
            indel_start_end = get_indel_start_end_pos(ll[2])
            indel_type = ll[2].split(":")[2]
            indel_pos_dic.setdefault(rname, {})[indel_start_end] = indel_type
    return indel_pos_dic

def get_nearest_indel(relative_pos, indel_dic):
    min_distace = 1e100
    min_direction = ""
    min_indel_type = ""

    for indel_start_end in indel_dic:
        indel_start = int(indel_start_end.split("\t")[0])
        indel_end = int(indel_start_end.split("\t")[1])
        indel_type = "\t".join(indel_dic[indel_start_end].split("_"))
        if relative_pos >= indel_start and relative_pos <= indel_end:
            min_direction = "inside"
            min_distace = min(abs(relative_pos - indel_start), abs(relative_pos - indel_end))
            min_indel_type = indel_type
            return min_direction, min_distace, min_indel_type

        elif relative_pos > indel_start and relative_pos > indel_end:
            direction = "right"
            distance = relative_pos - indel_end
        elif relative_pos < indel_start and relative_pos < indel_end:
            direction = "left"
            distance = indel_start - relative_pos
        
        if distance <= min_distace:
            min_distace = distance
            min_direction = direction
            min_indel_type = indel_type
    return min_direction, min_distace, min_indel_type


def read_snv(read_specific_snv, indel_pos_dic, out):
    with open(read_specific_snv, "r") as inf, open(out, "w") as final:
        for line in inf:
            ll = line.rstrip("\n").split("\t")
            rname = ll[0].split(":")[0]
            relative_pos = int(ll[1])

            if rname in indel_pos_dic:
                indel_dic = indel_pos_dic[rname]
                min_direction, min_distace, min_indel_type = get_nearest_indel(relative_pos, indel_dic)
                
                out_line = "\t".join(ll[0].split(":")) + "\t" + ll[1] + "\t" + ll[2] + "\t" + ll[3] + "\t" + str(min_direction) + "\t" + str(min_distace) + "\t" + str(min_indel_type)
                final.write(out_line + "\n")

if __name__ == "__main__":
    main()
