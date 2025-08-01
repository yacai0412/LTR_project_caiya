import sys

def main():
    inreadpos = sys.argv[1]
    inbed = sys.argv[2]

    out0 = inreadpos + "0"
    out1 = inreadpos + "1"

    sep_seq_start_0_base_dic = read_in_bed(inbed)
    read_inreadpos_convert_sep_pos(inreadpos, sep_seq_start_0_base_dic, out0, out1)
    

def read_in_bed(inbed):
    sep_seq_start_0_base_dic = {}
    with open(inbed, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            sep_rname = ll[3]
            start = int(ll[1])
            sep_seq_start_0_base_dic[sep_rname] = start
    return sep_seq_start_0_base_dic

def read_inreadpos_convert_sep_pos(inreadpos, sep_seq_start_0_base_dic, out0, out1):
    with open(inreadpos, "r") as inf, open(out0, "w") as final0, open(out1, "w") as final1:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")
            sep_rname = ll[0]
            sep_pos = int(ll[1])
            if sep_rname in sep_seq_start_0_base_dic:
                sep_start = sep_seq_start_0_base_dic[sep_rname]
                whole_read_name = "_".join(sep_rname.split("_")[:-1])
                whole_read_pos = sep_start + sep_pos
                outline1 = whole_read_name + "\t" + str(whole_read_pos)
                final1.write(outline1 + "\n")

                outline0 = line + "\t" + outline1
                final0.write(outline0 + "\n")
            else:
                print(sep_rname + " read not found in the bed" )


if __name__ == "__main__":
    main()
