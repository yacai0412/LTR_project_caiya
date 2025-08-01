import sys


def main():
    in_relativepos = sys.argv[1] # AAC.snv.1bp.inside_indel.refpos.relativepos

    out1 = in_relativepos + "1"

    realtive_pos_count_dic, total_count = read_relative_pos_file(in_relativepos)
    export_relative_pos_out1(realtive_pos_count_dic, total_count, out1)


def create_relativepos_dic():
    realtive_pos_count_dic = {}
    for j in ["DEL", "INS", "indel"]:
        realtive_pos_count_dic.setdefault(0, {})[j] = 0
    
        for i in range(10000):
            realtive_pos_count_dic.setdefault(i + 1, {})[j] = 0
            realtive_pos_count_dic.setdefault(-1 * (i + 1), {})[j] = 0

    return realtive_pos_count_dic

def read_relative_pos_file(in_relativepos):
    realtive_pos_count_dic = create_relativepos_dic()
    total_count = 0
    with open(in_relativepos, "r") as inf:
        for line in inf:
            total_count += 1

            line = line.rstrip("\n")
            ll = line.split("\t")
            relative_pos = int(ll[-2])
            indel_type = ll[-1]

            if relative_pos in realtive_pos_count_dic and indel_type in realtive_pos_count_dic[relative_pos]:
                realtive_pos_count_dic[relative_pos][indel_type] += 1
                realtive_pos_count_dic[relative_pos]["indel"] += 1
    return realtive_pos_count_dic, total_count

def export_relative_pos_out1(realtive_pos_count_dic, total_count, out1):
    with open(out1, "w") as final1:
        for j in ["DEL", "INS", "indel"]:
            outline0 = "0\t0\t" + j + "\t" + str(realtive_pos_count_dic[0][j]) + "\t" + str(realtive_pos_count_dic[0][j] / total_count)
            final1.write(outline0 + "\n")

            for i in range(10000):
                outline_right = "right\t" + str(i + 1) + "\t" + j + "\t" + str(realtive_pos_count_dic[i + 1][j]) + "\t" + str(realtive_pos_count_dic[i + 1][j] / total_count)
                outline_left = "left\t" + str(-1 * (i + 1)) + "\t" + j + "\t" + str(realtive_pos_count_dic[-1 * (i + 1)][j]) + "\t" + str(realtive_pos_count_dic[-1 * (i + 1)][j] / total_count)

                final1.write(outline_right + "\n")
                final1.write(outline_left + "\n")


if __name__ == "__main__":
    main()
